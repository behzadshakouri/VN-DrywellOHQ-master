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

#ifdef Behzad
    const string path="/home/behzad/Projects/VN Drywell_Models/";
    const string ohq_r="/home/behzad/Projects/OpenHydroQual/resources/";
#elif PowerEdge
    const string path="/mnt/3rd900/Projects/VN Drywell_Models/";
    const string ohq_r="/mnt/3rd900/Projects/OpenHydroQual/resources/";
#elif Arash
    const string path="/home/arash/Projects/VN Drywell_Models/";
    const string ohq_r="/home/arash/Projects/OpenHydroQual/resources/";
#elif SligoCreek
    const string path="/media/arash/E/Projects/VN Drywell_Models/";
    const string ohq_r="/media/arash/E/Projects/OpenHydroQual/resources/";
#endif

using namespace std;

int main(int argc, char *argv[])
{
    // ============================================================
    //   SIMULATION CONFIGURATION
    // ============================================================

    bool Model_Creator = 1; // 1 for using modelcreator, and 0 for loading saved Json file; Set it for every simulation

    SimulationConfig simcfg;
    RainConfig raincfg;

    // (1) Which rainfall dataset?
    // 1–3 = LA datasets, 4 = Pacoima, 5 = Synthetic
    raincfg.rain_data = 5; // Set to Pacoima

    // (2) Simulation batch index
    int Simulation_num = 1;          // Set for each run
    double Simulation_days = 2;    // Each run window

    // (3) Base time
    double Base_start;
    if (raincfg.rain_data < 4)
    Base_start = 43750;      // LA base
    else if (raincfg.rain_data == 4)
    Base_start = 43466;       // Pacoima base
    else if (raincfg.rain_data == 5)
    Base_start = 45763;       // Synthetic base

    double Base_end   = Base_start + Simulation_days;

    simcfg.Base_start = Base_start;
    simcfg.Base_end = Base_end;

    if (Model_Creator) {
        simcfg.start_time = simcfg.Base_start;
        simcfg.end_time   = simcfg.Base_end;
    } else {
        double shift = Simulation_days * (Simulation_num - 1);
        simcfg.start_time = simcfg.Base_start + shift;
        simcfg.end_time   = simcfg.Base_end   + shift;
    }

    // Solver limits (can override from main)
    simcfg.maximum_time_allowed = 10 * 86400;              // 10 days
    simcfg.maximum_matrix_inversions = 10 * 200000;        // 10× default

    cout << "=== Simulation Window ===\n";
    cout << "Start = " << simcfg.start_time << "\n";
    cout << "End   = " << simcfg.end_time << "\n";
    cout << "Model_Creator = " << Model_Creator << "\n\n";

    // ============================================================
    //   FIELD GENERATOR (soil stochasticity)
    // ============================================================

    FieldGenerator gen(200, 42); // Grid number and seed
    gen.setPDFMode(FieldGenerator::pdfmode::parametric); // parametric is new pdfmode

    // Set grid spacing to 0.5 meters
    gen.setDx(0.5);

    // Get current grid spacing
    double currentDx = gen.getDx();
    std::cout << "Grid spacing: " << currentDx << " m\n";

    // Input CSV and output folder
    std::string csvFile   = std::string(path) + "Soil retention params vs depth.csv";
    std::string outPrefix = std::string(path);

    std::vector<double> testDistances = {};
    std::vector<std::string> parameters = {"alpha", "n", "theta_s", "theta_r", "Ksat"};

    for (const auto &param : parameters) {
        generateAndAnalyzeField(gen, csvFile, param, 1, testDistances, outPrefix);
    }

    // ============================================================
    //   BUILD OR LOAD MODEL
    // ============================================================

    model_parameters mp;

    System *system = new System();
    system->SetWorkingFolder(path);

    ModelCreator ModCreate;

    if (Model_Creator)
    {
        cout << "Creating model...\n";
        ModCreate.Create(mp, system, &gen, raincfg, simcfg);
        cout << "Model build complete.\n\n";
    }
    else
    {
        cout << "Loading model...\n";
        system->ReadSystemSettingsTemplate(ohq_r + "settings.json");
        system->LoadfromJson(QString::fromStdString(path + "Model.json"));
        system->AddSolveVariableOrder("Storage");

        system->SetSettingsParameter("simulation_start_time", simcfg.start_time);
        system->SetSettingsParameter("simulation_end_time", simcfg.end_time);
        system->SetSystemSettings();

        cout << "Model loaded.\n\n";
    }

    // ============================================================
    //   SAVE SCRIPT + INITIALIZE
    // ============================================================

    system->SetSilent(false);
    cout<<"Saving"<<endl;
    system->SavetoScriptFile(path + "/CreatedModel.ohq");

    system->SavetoJson(path + "Model_test.json",system->addedtemplates, false, true );

    cout<<"Solving ..."<<endl;
    system->CalcAllInitialValues();

    // ============================================================
    //   SOIL PARAMETER SNAPSHOTS (VTP)
    // ============================================================

    ResultGrid resgrid_Ksat("K_sat_original",system);
    resgrid_Ksat.WriteToVTPSnapShot("Ksat",path + "Ksat.vtp",0,0);

    ResultGrid resgrid_alpha("alpha",system);
    resgrid_alpha.WriteToVTPSnapShot("alpha",path + "alpha.vtp",0,0);

    ResultGrid resgrid_n("n",system);
    resgrid_n.WriteToVTPSnapShot("n",path + "n.vtp",0,0);

    ResultGrid resgrid_theta_s("theta_sat",system);
    resgrid_theta_s.WriteToVTPSnapShot("theta_s",path + "theta_s.vtp",0,0);

    ResultGrid resgrid_theta_r("theta_res",system);
    resgrid_theta_r.WriteToVTPSnapShot("theta_r",path + "theta_r.vtp",0,0);

    vector<string> quantities = {"K_sat_original","alpha", "n","theta_sat", "theta_res"};
    double dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    ResultGrid::Make3DVTK(quantities,dr,system,path + "3D_model.vtm");

    // ============================================================
    //   SOLVE
    // ============================================================

    auto start = std::chrono::high_resolution_clock::now();
    system->Solve(false, false);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Solve took " << elapsed.count()/3600 << " hrs\n";

    // ============================================================
    //   SAVE MODEL + OUTPUTS
    // ============================================================

    system->SavetoJson(path + "Model.json",system->addedtemplates, false, true );
    system->SavetoJson(path + "Model_check.json",system->addedtemplates, true, true );
    cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() <<"'"<<endl;

    TimeSeriesSet<double> uniformoutput_LR = system->GetOutputs().make_uniform(1);
    TimeSeriesSet<double> uniformoutput_HR = system->GetOutputs().make_uniform(0.1);
    uniformoutput_LR.write(system->GetWorkingFolder() + "Output_LR.txt");
    uniformoutput_HR.write(system->GetWorkingFolder() + system->OutputFileName());
    cout<<"Getting results into grid"<<endl;

    // ============================================================
    //   RESULT GRIDS (Moisture + Age)
    // ============================================================

    double start_counter;
    if (Model_Creator)
        start_counter = 0;
    else
        start_counter = simcfg.start_time - simcfg.Base_start;

    ResultGrid resgrid(uniformoutput_LR,"theta",system);
    cout<<"Writing VTPs"<<endl;
    resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"Moisture/"+"moisture.vtp",0,start_counter);
    resgrid.write(system->GetWorkingFolder()+"theta_results.csv");

    ResultGrid resgrid_age(uniformoutput_LR,"meanagetracer:concentration",system);
    cout<<"Writing age VTPs"<<endl;
    resgrid_age.WriteToVTP("Mean Age",system->GetWorkingFolder()+"Moisture/"+"mean_age.vtp",0,start_counter);
    resgrid_age.write(system->GetWorkingFolder()+"age_results.csv");

    //  =========================================================
    //  ERT plots / snapshots @ measurement times
    //  =========================================================

    bool Use_ResultGrid_ERT = true;

    int nz_uw_n = (raincfg.rain_data == 5) ? mp.nz_uw_n : mp.nz_uw;

    // Only ERT-3 and ERT-5 for now (others commented out)
    std::vector<BoreholeSpec> boreholes = {
        {"ERT-3",  6.7979,  system->GetWorkingFolder() + "obs/ERT-3_obs.csv", 0.0},
        {"ERT-5",  4.3974,  system->GetWorkingFolder() + "obs/ERT-5_obs.csv", 0.0},

        // {"ERT-2",  0.9381,  system->GetWorkingFolder() + "obs/ERT-2_obs.csv", 0.0},
        // {"ERT-4",  4.1679,  system->GetWorkingFolder() + "obs/ERT-4_obs.csv", 0.0},
        // {"ERT-1", 10.0,     system->GetWorkingFolder() + "obs/ERT-1_obs.csv", 0.0},
    };

    BoreholeExportOptions opt;
    opt.nz_uw_n          = nz_uw_n;
    opt.out_dir          = system->GetWorkingFolder();
    opt.model_theta_var  = "theta";
    opt.allow_missing_obs = true;
    opt.use_resultgrid   = Use_ResultGrid_ERT;

    // Your ERT measurement times (Excel serial days; Base_start=45763 => 2025-04-16 00:00)
    const std::vector<double> ert_times = {
        45763.588194444, // 4/16/2025 14:07
        45763.682638889, // 4/16/2025 16:23
        45763.770833333, // 4/16/2025 18:30
        45764.673611111, // 4/17/2025 16:10
        45764.729166667  // 4/17/2025 17:30
    };

    // File-friendly token (replace '.' with 'p')
    auto time_token = [&](double t) -> std::string
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << t;   // stable filenames
        std::string s = oss.str();
        std::replace(s.begin(), s.end(), '.', 'p');
        return s;
    };

    // Find nearest time index in uniformoutput_LR to a target serial-day time
    auto nearest_time_index = [&](double target_t) -> unsigned int
    {
        if (uniformoutput_LR.size() == 0) return 0;

        // Use first series as reference time axis
        const auto &ts0 = uniformoutput_LR[0];
        unsigned int n = (unsigned int)ts0.size();
        if (n == 0) return 0;

        unsigned int best_i = 0;
        double best_dt = std::numeric_limits<double>::infinity();

        for (unsigned int i = 0; i < n; ++i)
        {
            double t = ts0.getTime(i);
            double dt = std::abs(t - target_t);
            if (dt < best_dt) { best_dt = dt; best_i = i; }
        }
        return best_i;
    };

    // Export snapshots/CSVs for each ERT measurement time
    for (double t_ert : ert_times)
    {
        unsigned int idx = nearest_time_index(t_ert);

        // IMPORTANT: prevent overwriting outputs by changing prefix per time
        // (ERT_comparison writes: out_dir + file_prefix + bh.name + file_suffix)
        std::string tag = "t" + time_token(t_ert);
        opt.file_prefix = "ERTsnap_" + tag + "_";   // e.g., ERTsnap_t45763p588194_...
        opt.file_suffix = ".csv";                  // keep default format
        opt.time_index  = idx;

        std::cout << "\n=== ERT snapshot export ===\n";
        std::cout << "Target ERT time = " << std::setprecision(12) << t_ert << "\n";
        std::cout << "Nearest model index = " << idx << "\n";
        std::cout << "Output prefix = " << opt.file_prefix << "\n";

        export_borehole_theta_comparison_csvs(
            uniformoutput_LR, mp, system,
            boreholes, opt
        );
    }


    // ============================================================
    //   WELL DEPTHS
    // ============================================================

    vector<string> well_block_c; well_block_c.push_back("Well_c");
    ResultGrid well_depth_c = ResultGrid(uniformoutput_HR,well_block_c,"depth");
    well_depth_c.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_c.csv");

    // ============================================================
    //   MEAN AGE AT GROUNDWATER + RECHARGE
    // ============================================================

    vector<string> age_tracer_locations;
    vector<string> gw_rechage_locations;
    for (unsigned int i = 0; i<=ModCreate.ModelParameters().nr_uw; i++)
    {
        age_tracer_locations.push_back("Soil-uw ("+ aquiutils::numbertostring(i) +"$"+ aquiutils::numbertostring(ModCreate.ModelParameters().nz_uw-1)+")");
        gw_rechage_locations.push_back("Soil to Groundwater (" + aquiutils::numbertostring(i) + ")");
    }
    ResultGrid age_tracer_res = ResultGrid(uniformoutput_HR,age_tracer_locations, "meanagetracer:concentration");
    ResultGrid flow_to_gw_res = ResultGrid(uniformoutput_HR,gw_rechage_locations, "flow");

    TimeSeries<double> mean_age;
    TimeSeries<double> groundwater_recharge;
    for (unsigned int j = 0; j<age_tracer_res.maxnumpoints(); j++)
    {
        double sumprod = 0;
        double sumflow = 0;
        for (unsigned int i = 0; i<ModCreate.ModelParameters().nr_uw; i++)
        {
            sumprod += age_tracer_res[i].getValue(j)*flow_to_gw_res[i].getValue(j);
            sumflow += flow_to_gw_res[i].getValue(j);
        }
        mean_age.append(age_tracer_res[0].getTime(j),sumprod/sumflow);
        groundwater_recharge.append(age_tracer_res[0].getTime(j),sumflow);
    }

    mean_age.writefile(system->GetWorkingFolder()+"mean_age_at_gw.csv");
    groundwater_recharge.writefile(system->GetWorkingFolder()+"gw_recharge.csv");

    vector<string> well_block_g; well_block_g.push_back("Well_g");
    ResultGrid well_depth_g = ResultGrid(uniformoutput_HR,well_block_g,"depth");
    well_depth_g.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_g.csv");

    // ============================================================
    //   SAVE STATE VARIABLES
    // ============================================================

    system->SaveStateVariableToJson("meanagetracer:concentration",path + "Age.json");
    system->SaveStateVariableToJson("theta",path + "theta.json");

    return 0;
}
