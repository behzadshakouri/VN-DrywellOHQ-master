#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"
#include "ostream"
#include "fieldgenerator.h"
#include "FieldGenHelper.h"
#include "ERT_comparison.h"

#include <chrono>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef Behzad
    const std::string path  = "/home/behzad/Projects/VN Drywell_Models/";
    const std::string ohq_r = "/home/behzad/Projects/OpenHydroQual/resources/";
#elif PowerEdge
    const std::string path  = "/mnt/3rd900/Projects/VN Drywell_Models/";
    const std::string ohq_r = "/mnt/3rd900/Projects/OpenHydroQual/resources/";
#elif Arash
    const std::string path  = "/home/arash/Projects/VN Drywell_Models/";
    const std::string ohq_r = "/home/arash/Projects/OpenHydroQual/resources/";
#elif SligoCreek
    const std::string path  = "/media/arash/E/Projects/VN Drywell_Models/";
    const std::string ohq_r = "/media/arash/E/Projects/OpenHydroQual/resources/";
#endif

using namespace std;

int main(int argc, char *argv[])
{
    // ============================================================
    //   SIMULATION CONFIGURATION
    // ============================================================

    bool Model_Creator = 1;   // 1 = build model and solve; 0 = load model and re-use saved outputs
    bool Run_Solve     = Model_Creator;

    SimulationConfig simcfg;
    RainConfig raincfg;

    // ------------------------------------------
    // ERT analysis run = synthetic rainfall (5)
    // ------------------------------------------
    raincfg.rain_data = 5;   // Synthetic rainfall

    int    Simulation_num  = 1;
    double Simulation_days = 2;

    // Base_start is an Excel-day reference for your rainfall datasets
    double Base_start;
    if (raincfg.rain_data < 4)
        Base_start = 43750;
    else if (raincfg.rain_data == 4)
        Base_start = 43466;
    else
        Base_start = 45763;

    double Base_end = Base_start + Simulation_days;

    simcfg.Base_start = Base_start;
    simcfg.Base_end   = Base_end;

    if (Model_Creator)
    {
        // New run: run from the base window
        simcfg.start_time = simcfg.Base_start;
        simcfg.end_time   = simcfg.Base_end;
    }
    else
    {
        // Re-using outputs: optionally shift the window by Simulation_num
        double shift = Simulation_days * (Simulation_num - 1);
        simcfg.start_time = simcfg.Base_start + shift;
        simcfg.end_time   = simcfg.Base_end   + shift;
    }

    simcfg.maximum_time_allowed      = 10 * 86400;
    simcfg.maximum_matrix_inversions = 10 * 200000;

    // ============================================================
    //   Ksat Scale Factor (CLI)
    //
    //   unsaturated_soil_revised_model.json must use:
    //       K_sat = K_sat_original * K_sat_scale_factor
    //
    //   CLI examples:
    //       --ksat-scale 10
    //       --ksat-scale=10
    //
    //   Optional per-group overrides:
    //       --ksat-scale-g 2
    //       --ksat-scale-uw 0.5
    //
    //   Behavior:
    //     - default applies to ALL soils
    //     - if g/uw override is provided (finite >0), it overrides only that group
    // ============================================================
    simcfg.KsatScaleFactor    = parse_ksat_scale(argc, argv, 1.0);
    simcfg.KsatScaleFactor_g  = parse_ksat_scale_g(argc, argv, 1.0);
    simcfg.KsatScaleFactor_uw = parse_ksat_scale_uw(argc, argv, 5.0);

    // ============================================================
    //   INITIAL THETA MODE (AUTO for ERT analysis runs)
    // ============================================================

    const bool is_ert_analysis_run = (raincfg.rain_data == 5);

    const bool user_set_init_theta = cli_has_init_theta(argc, argv);
    InitThetaMode initThetaMode = parse_init_theta_mode(argc, argv);

    if (is_ert_analysis_run && !user_set_init_theta)
        initThetaMode = InitThetaMode::ERT_IDW_R;

    const bool useERTInit = (initThetaMode != InitThetaMode::Default);

    int mc_mode  = 0;
    int mc_which = 0;

    if (initThetaMode == InitThetaMode::ERT3_Only) { mc_which = 3; mc_mode = 0; }
    if (initThetaMode == InitThetaMode::ERT5_Only) { mc_which = 5; mc_mode = 0; }
    if (initThetaMode == InitThetaMode::ERT_IDW_R)  { mc_which = 0; mc_mode = 0; }
    if (initThetaMode == InitThetaMode::ERT_R_Avg)  { mc_which = 0; mc_mode = 1; }

    cout << "=== Simulation Window ===\n";
    cout << "Start = " << simcfg.start_time << "\n";
    cout << "End   = " << simcfg.end_time   << "\n";
    cout << "Model_Creator = " << Model_Creator << "\n";
    cout << "Run_Solve     = " << Run_Solve     << "\n";
    cout << "Rain data     = " << raincfg.rain_data
         << (is_ert_analysis_run ? " (ERT analysis)\n" : "\n");

    cout << "InitThetaMode = " << init_theta_mode_name(initThetaMode)
         << " (user_set=" << (user_set_init_theta ? 1 : 0)
         << ", UseERTInitialTheta=" << (useERTInit ? 1 : 0)
         << ", mc_mode=" << mc_mode
         << ", mc_which=" << mc_which << ")\n";

    cout << "KsatScaleFactor(default) = " << simcfg.KsatScaleFactor << "\n";
    cout << "KsatScaleFactor_g        = " << simcfg.KsatScaleFactor_g  << "\n";
    cout << "KsatScaleFactor_uw       = " << simcfg.KsatScaleFactor_uw << "\n\n";

    // ============================================================
    //   FIELD GENERATOR
    // ============================================================

    FieldGenerator gen(200, 42);
    gen.setPDFMode(FieldGenerator::pdfmode::parametric);
    gen.setDx(0.5);

    cout << "Grid spacing: " << gen.getDx() << " m\n";

    std::string csvFile   = path + "Soil retention params vs depth.csv";
    std::string outPrefix = path;

    std::vector<double> testDistances = {};
    std::vector<std::string> parameters = {"alpha", "n", "theta_s", "theta_r", "Ksat"};

    for (const auto &param : parameters)
        generateAndAnalyzeField(gen, csvFile, param, 1, testDistances, outPrefix);

    // ============================================================
    //   BUILD OR LOAD MODEL
    // ============================================================

    model_parameters mp;

    System *system = new System();
    system->SetWorkingFolder(path);

    ModelCreator ModCreate;

    ModCreate.UseERTInitialTheta   = useERTInit;
    ModCreate.ERTInitialThetaMode  = mc_mode;
    ModCreate.ERTInitialThetaWhich = mc_which;

    if (Model_Creator)
    {
        cout << "Creating model...\n";
        ModCreate.Create(mp, system, &gen, raincfg, simcfg);
        mp = ModCreate.ModelParameters();
        cout << "Model build complete.\n\n";
    }
    else
    {
        cout << "Loading model...\n";
        system->ReadSystemSettingsTemplate(ohq_r + "settings.json");
        system->LoadfromJson(QString::fromStdString(path + "Model.json"));
        system->AddSolveVariableOrder("Storage");

        system->SetSettingsParameter("simulation_start_time", simcfg.start_time);
        system->SetSettingsParameter("simulation_end_time",   simcfg.end_time);
        system->SetSystemSettings();

        mp = ModCreate.ModelParameters();
        cout << "Model loaded.\n\n";
    }

    // ============================================================
    //   SAVE SCRIPT + INITIALIZE
    // ============================================================

    system->SetSilent(false);

    cout << "Saving model files...\n";
    system->SavetoScriptFile(path + "/CreatedModel.ohq");
    system->SavetoJson(path + "Model_test.json", system->addedtemplates, false, true);

    cout << "CalcAllInitialValues ...\n";
    system->CalcAllInitialValues();

    // ============================================================
    //   SOLVE (optional)
    // ============================================================

    if (Run_Solve)
    {
        auto start = std::chrono::high_resolution_clock::now();
        system->Solve(false, false);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Solve took " << elapsed.count()/3600 << " hrs\n";
    }
    else
    {
        std::cout << "Run_Solve=0 -> skipping Solve(). Will read saved outputs from disk.\n";
    }

    // ============================================================
    //   LOAD OUTPUTS (if not solved) OR GET OUTPUTS (if solved)
    // ============================================================

    TimeSeriesSet<double> rawOutputs;

    if (Run_Solve)
    {
        rawOutputs = system->GetOutputs();
    }
    else
    {
        std::string outFile = system->GetWorkingFolder() + system->OutputFileName();
        std::string outLR   = system->GetWorkingFolder() + "Output_LR.txt";

        QString qOutFile = QString::fromStdString(outFile);
        QString qOutLR   = QString::fromStdString(outLR);

        if (QFileInfo::exists(qOutFile))
        {
            std::cout << "Reading saved outputs: " << outFile << "\n";
            rawOutputs.read(outFile);
        }
        else if (QFileInfo::exists(qOutLR))
        {
            std::cout << "Reading saved outputs: " << outLR << "\n";
            rawOutputs.read(outLR);
        }
        else
        {
            std::cerr << "ERROR: No saved output file found.\n";
            std::cerr << "Looked for:\n  " << outFile << "\n  " << outLR << "\n";
            return 1;
        }
    }

    // Uniform time grids
    TimeSeriesSet<double> uniformoutput_LR  = rawOutputs.make_uniform(1);
    TimeSeriesSet<double> uniformoutput_HR  = rawOutputs.make_uniform(0.1);

    // ERT grid: 1 minute in Excel-days
    double dt_ert_days = 1.0 / 1440.0;
    TimeSeriesSet<double> uniformoutput_ERT = rawOutputs.make_uniform(dt_ert_days);

    // ============================================================
    //   RESULT GRIDS
    // ============================================================

    double start_counter = Model_Creator ? 0 : (simcfg.start_time - simcfg.Base_start);

    ResultGrid resgrid(uniformoutput_LR, "theta", system);
    cout << "Writing VTPs\n";
    resgrid.WriteToVTP("Moisture_content",
                       system->GetWorkingFolder() + "Moisture/" + "moisture.vtp",
                       0, start_counter);
    resgrid.write(system->GetWorkingFolder() + "theta_results.csv");

    ResultGrid resgrid_age(uniformoutput_LR, "meanagetracer:concentration", system);
    cout << "Writing age VTPs\n";
    resgrid_age.WriteToVTP("Mean Age",
                           system->GetWorkingFolder() + "Moisture/" + "mean_age.vtp",
                           0, start_counter);
    resgrid_age.write(system->GetWorkingFolder() + "age_results.csv");

    // ============================================================
    //   ERT plots / snapshots @ measurement times
    // ============================================================
    export_ert_snapshots_at_measurement_times(
        uniformoutput_ERT, mp, system, raincfg,
        /*use_resultgrid_ert=*/true,
        /*allow_missing_obs=*/true
    );

    // ============================================================
    //   SAVE STATE VARIABLES
    // ============================================================

    system->SaveStateVariableToJson("meanagetracer:concentration", path + "Age.json");
    system->SaveStateVariableToJson("theta",                      path + "theta.json");

    return 0;
}
