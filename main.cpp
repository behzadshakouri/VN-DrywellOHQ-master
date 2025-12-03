#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"
#include "ostream"
#include "fieldgenerator.h"
#include "FieldGenHelper.h"

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

    bool Model_Creator = 1; // 1 for using modelcreator, and 0 for loading saved Json file; Set it for every simulation

    double Simulation_num = 8; // Simulation Number; Set it for every simulation
    double Simulation_days = 180; // Simulation Days for Json file [Constant]

    double Simulation_start_time_0=43466; // Simulation Start Date for Model_Creator [Pacoima station]
    double Simulation_end_time_0=45292; // Simulation End Date for Model_Creator [Pacoima station]
    //double Simulation_start_time_0 = 43750; // Simulation Start Date for Model_Creator
    //double Simulation_end_time_0 = 43930; // Simulation Start Date for Model_Creator
    //double Simulation_start_time_0 = 44864; // Simulation Start Date for Model_Creator [Test]
    //double Simulation_end_time_0 = 44866; // Simulation End Date for Model_Creator [Test]

    double Simulation_days_0 = Simulation_end_time_0 - Simulation_start_time_0; // Simulation Days for Model_Creator

    double Simulation_start_time; // Simulation Start Date
    double Simulation_end_time; // Simulation End Date

    // For Json file

    if (Model_Creator)
    {
        Simulation_start_time = Simulation_start_time_0; // Simulation Start Date
        Simulation_end_time = Simulation_end_time_0; // Simulation End Date
     }
    else
    {
        Simulation_start_time = Simulation_start_time_0 + Simulation_days * (Simulation_num - 1); // Simulation Start Date
        Simulation_end_time = Simulation_end_time_0 + Simulation_days * (Simulation_num - 1); // Simulation End Date
    }

    //omp_set_nested(0);          // Disable nested parallelism
    //omp_set_dynamic(0);         // Optional: disable dynamic thread adjustment

    // === Field Generator Test ===

    FieldGenerator gen(200, 42); // Grid number and seed
    gen.setPDFMode(FieldGenerator::pdfmode::parametric);

    // Set grid spacing to 0.5 meters
    gen.setDx(0.5);

    // Get current grid spacing
    double currentDx = gen.getDx();
    std::cout << "Grid spacing: " << currentDx << " m\n";

    // Input CSV and output folder
    std::string csvFile   = std::string(path) + "Soil retention params vs depth.csv";
    std::string outPrefix = std::string(path);

    // Distances to check auto-correlation
    //std::vector<double> testDistances = {0.0, 0.2, 0.5, 1.0, 2.0, 2.5, 5.0, 7.5, 10.0, 20.0, 50.0, 100.0};
    std::vector<double> testDistances = {};
    // Parameters to generate
    std::vector<std::string> parameters = {"alpha", "n", "theta_s", "theta_r", "Ksat"};

    // Loop through parameters
    for (const auto &param : parameters) {
        generateAndAnalyzeField(gen, csvFile, param, 1, testDistances, outPrefix); // enter correlation length
    }

    model_parameters mp;

    System *system=new System();

    ModelCreator ModCreate;

    system->SetWorkingFolder(path); // Should be modified according to the users directory

    if (Model_Creator)
    {
        cout<<"Creating model ..." <<endl;
        //ModCreate.Create(mp,system, &gen);
        ModCreate.Create(mp,system);
        cout<<"Creating model done..." <<endl;
    }
    else
    {
        system->ReadSystemSettingsTemplate(ohq_r+"settings.json");
        system->LoadfromJson(QString::fromStdString(path + "Model.json"));

        system->AddSolveVariableOrder("Storage"); // Flow solve

    cout<<"Model loaded..." <<endl;
    // Solve properties
        system->SetSettingsParameter("simulation_start_time",Simulation_start_time);
        system->SetSettingsParameter("simulation_end_time",Simulation_end_time);
        system->SetSystemSettings();
    }

    system->SetSilent(false);
    cout<<"Saving"<<endl;
    system->SavetoScriptFile(path + "/CreatedModel.ohq"); // Should be modified according to the users directory

    cout<<"Solving ..."<<endl;
    system->CalcAllInitialValues();

    // Soil params VTP
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


    // Solve
    auto start = std::chrono::high_resolution_clock::now();
    system->Solve(false, false);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Solve took " << elapsed.count()/3600 << " hrs\n";


    system->SavetoJson(path + "Model.json",system->addedtemplates, false, true );
    system->SavetoJson(path + "Model_check.json",system->addedtemplates, true, true );
    cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() +"'"<<endl;

    TimeSeriesSet<double> uniformoutput_LR = system->GetOutputs().make_uniform(1);
    TimeSeriesSet<double> uniformoutput_HR = system->GetOutputs().make_uniform(0.1);
    uniformoutput_LR.write(system->GetWorkingFolder() + "Output_LR.txt");
    uniformoutput_HR.write(system->GetWorkingFolder() + system->OutputFileName());
    cout<<"Getting results into grid"<<endl;

    double start_counter;
    if (Model_Creator)
        start_counter = 0;
    else
        start_counter = Simulation_start_time - Simulation_start_time_0;

    ResultGrid resgrid(uniformoutput_LR,"theta",system);
    cout<<"Writing VTPs"<<endl;
    resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"Moisture/"+"moisture.vtp",0,start_counter);
    resgrid.write(system->GetWorkingFolder()+"theta_results.csv");

    ResultGrid resgrid_age(uniformoutput_LR,"meanagetracer:concentration",system);
    cout<<"Writing age VTPs"<<endl;
    resgrid_age.WriteToVTP("Mean Age",system->GetWorkingFolder()+"Moisture/"+"mean_age.vtp",0,start_counter);
    resgrid_age.write(system->GetWorkingFolder()+"age_results.csv");

    vector<string> well_block_c; well_block_c.push_back("Well_c");
    ResultGrid well_depth_c = ResultGrid(uniformoutput_HR,well_block_c,"depth");
    well_depth_c.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_c.csv");

    vector<string> age_tracer_locations;
    vector<string> gw_rechage_locations;
    for (unsigned int i = 0; i<=ModCreate.ModelParameters().nr_uw; i++)
    {   age_tracer_locations.push_back("Soil-uw ("+ aquiutils::numbertostring(i) +"$"+ aquiutils::numbertostring(ModCreate.ModelParameters().nz_uw-1)+")");
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

    system->SaveStateVariableToJson("meanagetracer:concentration",path + "Age.json");
    system->SaveStateVariableToJson("theta",path + "theta.json");

    return 0;

}
