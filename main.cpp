#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"
#include "ostream"
#include "fieldgenerator.h"

#ifdef PowerEdge
#define PATH "/mnt/3rd900/Projects/VN Drywell_Models/"
#elif Arash
#define PATH "/home/arash/Projects/VN Drywell_Models/"
#endif

using namespace std;

// Helper function to generate, analyze, and save a field
void generateAndAnalyzeField(FieldGenerator &gen,
                             const std::string &csvPath,
                             const std::string &paramName,
                             double correlationLength,
                             const std::vector<double> &testDistances,
                             const std::string &outputPrefix)
{
    std::string fieldName = paramName + "_normal";

    std::cout << "Reading data from '" << csvPath << "'" << std::endl;

    // Generate field
    gen.generateFieldFromMeasurementData(csvPath, paramName, fieldName, correlationLength);

    // Save CDFs
    gen.getMeasuredCDFs().write(outputPrefix + "Params_CDFs.txt");

    // Print stats
    std::cout << "Mean: " << gen.mean(fieldName) << std::endl;
    std::cout << "Standard Deviation: " << gen.standardDeviation(fieldName) << std::endl;
    std::cout << "Field size: " << 1000 << " points" << std::endl;
    std::cout << "Grid spacing: " << gen.getDx() << " m" << std::endl;
    std::cout << "Total length: " << 1000 * gen.getDx() << " m" << std::endl;

    // Auto-correlation analysis
    std::cout << "\n=== AUTO-CORRELATION ANALYSIS for " << paramName << " ===" << std::endl;
    std::cout << "Distance (m)  |  Correlation  |  Theoretical*" << std::endl;
    std::cout << "------------- | ------------- | -------------" << std::endl;

    for (double dist : testDistances) {
        double empirical = gen.autoCorrelation(fieldName, dist);
        double theoretical = std::exp(-dist / correlationLength);
        std::cout << std::setw(10) << dist << "   |  "
                  << std::setw(10) << empirical << "   |  "
                  << std::setw(10) << theoretical << std::endl;
    }

    // Save fields
    std::string normalFile = outputPrefix + fieldName + "_field.csv";
    std::string realFile   = outputPrefix + paramName + "_field.csv";

    if (gen.saveFieldToCSV(fieldName, normalFile))
        std::cout << "Saved " << normalFile << " successfully!\n";
    else
        std::cout << "Failed to save " << normalFile << "\n";

    if (gen.saveFieldToCSV(paramName, realFile))
        std::cout << "Saved " << realFile << " successfully!\n";
    else
        std::cout << "Failed to save " << realFile << "\n";
}

int main(int argc, char *argv[])
{

    // === Field Generator Test ===
    FieldGenerator gen(100, 42);

    // Set grid spacing to 0.5 meters
    gen.setDx(0.5);

    // Get current grid spacing
    double currentDx = gen.getDx();
    std::cout << "Grid spacing: " << currentDx << " m\n";

    // Input CSV and output folder
    std::string csvFile   = std::string(PATH) + "Soil retention params vs depth.csv";
    std::string outPrefix = std::string(PATH);

    // Distances to check auto-correlation
    std::vector<double> testDistances = {0.0, 0.2, 0.5, 1.0, 2.0, 2.5, 5.0, 7.5, 10.0, 20.0};

    // === Generate fields for all parameters ===
    generateAndAnalyzeField(gen, csvFile, "Ksat",    20, testDistances, outPrefix);
    generateAndAnalyzeField(gen, csvFile, "alpha",   20, testDistances, outPrefix);
    generateAndAnalyzeField(gen, csvFile, "n",       20, testDistances, outPrefix);
    generateAndAnalyzeField(gen, csvFile, "theta_s", 20, testDistances, outPrefix);
    generateAndAnalyzeField(gen, csvFile, "theta_r", 20, testDistances, outPrefix);

    return 0;

/*
    // === Ksat generation ===
    std::cout<<"Reading data from '" <<std::string(PATH) + "Soil retention params vs depth.csv'"<<std::endl;
        gen.generateFieldFromMeasurementData(
        std::string(PATH) + "Soil retention params vs depth.csv",
        "Ksat",                    // Parameter name in the CSV
        "Ksat_normal",           // Field name to generate
        20                       // Correlation length
        );

    gen.getMeasuredCDFs().write(std::string(PATH) + "Params_CDFs.txt");

    std::cout << "Mean: " << gen.mean("Ksat_normal") << std::endl;
    std::cout << "Standard Deviation: " << gen.standardDeviation("Ksat_normal") << std::endl;
    std::cout << "Field size: " << 1000 << " points" << std::endl;
    std::cout << "Grid spacing: " << gen.getDx() << " m" << std::endl;
    std::cout << "Total length: " << 1000 * gen.getDx() << " m" << std::endl;

    // === AUTO-CORRELATION ANALYSIS ===
    std::cout << "\n=== AUTO-CORRELATION ANALYSIS ===" << std::endl;

    // Test correlation at various distances
    std::vector<double> testDistances = {0.0, 0.2, 0.5, 1.0, 2.0, 2.5, 5.0, 7.5, 10.0, 20.0};

    std::cout << "Distance (m)  |  Correlation  |  Theoretical*" << std::endl;
    std::cout << "------------- | ------------- | -------------" << std::endl;

    for (double dist : testDistances) {
        double empirical = gen.autoCorrelation("Ksat_normal", dist);
        double theoretical = std::exp(-dist / 20.0);  // Theoretical exponential decay

        std::cout << std::setw(10) << dist << "   |  "
                  << std::setw(10) << empirical << "   |  "
                  << std::setw(10) << theoretical << std::endl;
    }

    if (gen.saveFieldToCSV("Ksat_normal", std::string(PATH) + "Ksat_normal_field.csv")) {
        std::cout << "Field saved successfully!\n";
    } else {
        std::cout << "Failed to save field\n";
    }

    if (gen.saveFieldToCSV("Ksat", std::string(PATH) + "Ksat_field.csv")) {
        std::cout << "Field saved successfully!\n";
    } else {
        std::cout << "Failed to save field\n";
    }

    // === alpha generation ===

    std::cout<<"Reading data from '" <<std::string(PATH) + "Soil retention params vs depth.csv'"<<std::endl;
    gen.generateFieldFromMeasurementData(
        std::string(PATH) + "Soil retention params vs depth.csv",
        "alpha",                    // Parameter name in the CSV
        "alpha_normal",           // Field name to generate
        20                       // Correlation length
        );

    gen.getMeasuredCDFs().write(std::string(PATH) + "Params_CDFs.txt");

    std::cout << "Mean: " << gen.mean("alpha_normal") << std::endl;
    std::cout << "Standard Deviation: " << gen.standardDeviation("alpha_normal") << std::endl;
    std::cout << "Field size: " << 1000 << " points" << std::endl;
    std::cout << "Grid spacing: " << gen.getDx() << " m" << std::endl;
    std::cout << "Total length: " << 1000 * gen.getDx() << " m" << std::endl;

    // === AUTO-CORRELATION ANALYSIS ===
    std::cout << "\n=== AUTO-CORRELATION ANALYSIS ===" << std::endl;

    // Test correlation at various distances
    std::vector<double> testDistances2 = {0.0, 0.2, 0.5, 1.0, 2.0, 2.5, 5.0, 7.5, 10.0, 20.0};

    std::cout << "Distance (m)  |  Correlation  |  Theoretical*" << std::endl;
    std::cout << "------------- | ------------- | -------------" << std::endl;

    for (double dist : testDistances2) {
        double empirical = gen.autoCorrelation("alpha_normal", dist);
        double theoretical = std::exp(-dist / 20.0);  // Theoretical exponential decay

        std::cout << std::setw(10) << dist << "   |  "
                  << std::setw(10) << empirical << "   |  "
                  << std::setw(10) << theoretical << std::endl;
    }

    if (gen.saveFieldToCSV("alpha_normal", std::string(PATH) + "alpha_normal_field.csv")) {
        std::cout << "Field saved successfully!\n";
    } else {
        std::cout << "Failed to save field\n";
    }

    if (gen.saveFieldToCSV("alpha", std::string(PATH) + "alpha_field.csv")) {
        std::cout << "Field saved successfully!\n";
    } else {
        std::cout << "Failed to save field\n";
    }

    return 0;
*/

    model_parameters mp;

    System *system=new System();
    ModelCreator ModCreate;
    cout<<"Creating model ..." <<endl;
    ModCreate.Create(mp,system, &gen);
    cout<<"Creating model done..." <<endl;

#ifdef PowerEdge
    string path = "/mnt/3rd900/Projects/VN Drywell_Models/";
#elif Arash
    string path = "/home/arash/Projects/VN Drywell_Models/";
#elif SligoCreek
    string path = "/media/arash/E/Projects/VN Drywell_Models/";
#endif

    system->SetWorkingFolder(path); // Should be modified according to the users directory
    system->SetSilent(false);
    cout<<"Saving"<<endl;
    system->SavetoScriptFile(path + "/CreatedModel.ohq"); // Should be modified according to the users directory

    cout<<"Solving ..."<<endl;
    system->CalcAllInitialValues();

    // Soil params VTP
    ResultGrid resgrid_Ksat("K_sat_original",system);
    resgrid_Ksat.WriteToVTP("Ksat",path + "Ksat.vtp",0);

    ResultGrid resgrid_alpha("alpha",system);
    resgrid_alpha.WriteToVTP("alpha",path + "alpha.vtp",0);

    ResultGrid resgrid_n("n",system);
    resgrid_n.WriteToVTP("n",path + "n.vtp",0);

    ResultGrid resgrid_theta_s("theta_sat",system);
    resgrid_theta_s.WriteToVTP("theta_s",path + "theta_s.vtp",0);

    ResultGrid resgrid_theta_r("theta_res",system);
    resgrid_theta_r.WriteToVTP("theta_r",path + "theta_r.vtp",0);

    vector<string> quantities = {"K_sat_original","alpha", "n","theta_sat", "theta_res"};
    double dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    ResultGrid::Make3DVTK(quantities,dr,system,path + "3D_model.vtm");


    // Solve
    system->Solve();
    system->SavetoJson(path + "/Model.json",system->addedtemplates, true, true );
    cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() +"'"<<endl;

    TimeSeriesSet<double> uniformoutput_LR = system->GetOutputs().make_uniform(1);
    TimeSeriesSet<double> uniformoutput_HR = system->GetOutputs().make_uniform(0.1);
    uniformoutput_HR.write(system->GetWorkingFolder() + system->OutputFileName());
    cout<<"Getting results into grid"<<endl;

    ResultGrid resgrid(uniformoutput_LR,"theta",system);
    cout<<"Writing VTPs"<<endl;
    resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"Moisture/"+"moisture.vtp");


    vector<string> well_block_c; well_block_c.push_back("Well_c");
    ResultGrid well_depth_c = ResultGrid(uniformoutput_HR,well_block_c,"depth");
    well_depth_c.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_c.csv");

    vector<string> well_block_g; well_block_g.push_back("Well_g");
    ResultGrid well_depth_g = ResultGrid(uniformoutput_HR,well_block_g,"depth");
    well_depth_g.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_g.csv");

    return 0;

}
