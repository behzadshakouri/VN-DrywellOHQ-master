#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"
#include "ostream"


#define PATH "/mnt/3rd900/Projects/VN Drywell_Models/"

using namespace std;

int main(int argc, char *argv[])
{

    model_parameters mp;

    System *system=new System();
    ModelCreator ModCreate;
    cout<<"Creating model ..." <<endl;
    ModCreate.Create(mp,system);
    cout<<"Creating model done..." <<endl;

#ifdef PowerEdge
    string path = "/mnt/3rd900/Projects/VN Drywell_Models/";
#elif Arash
    string path = "/home/arash/Projects/VN Drywell_Models";
#endif

    system->SetWorkingFolder(path); // Should be modified according to the users directory
    system->SetSilent(false);
    cout<<"Saving"<<endl;
    system->SavetoScriptFile(path + "/CreatedModel.ohq"); // Should be modified according to the users directory

    cout<<"Solving ..."<<endl;
    system->CalcAllInitialValues();

    // Soil params VTP
    ResultGrid resgrid_Ksat("K_sat_original",system);
    resgrid_Ksat.WriteToVTP("Ksat",path + "/Ksat.vtp",0);

    ResultGrid resgrid_alpha("alpha",system);
    resgrid_alpha.WriteToVTP("alpha",path + "/alpha.vtp",0);

    ResultGrid resgrid_n("n",system);
    resgrid_n.WriteToVTP("n",path + "/n.vtp",0);

    ResultGrid resgrid_theta_s("theta_sat",system);
    resgrid_theta_s.WriteToVTP("theta_s",path + "/theta_s.vtp",0);

    ResultGrid resgrid_theta_r("theta_res",system);
    resgrid_theta_r.WriteToVTP("theta_r",path + "/theta_r.vtp",0);

    vector<string> quantities = {"K_sat_original","alpha", "n","theta_sat", "theta_res"};
    double dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    ResultGrid::Make3DVTK(quantities,dr,system,path + "/3D_model.vtm");


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
    resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"/Moisture/"+"moisture.vtp");


    vector<string> well_block_c; well_block_c.push_back("Well_c");
    ResultGrid well_depth_c = ResultGrid(uniformoutput_HR,well_block_c,"depth");
    well_depth_c.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_c.csv");

    vector<string> well_block_g; well_block_g.push_back("Well_g");
    ResultGrid well_depth_g = ResultGrid(uniformoutput_HR,well_block_g,"depth");
    well_depth_g.Sum().writefile(system->GetWorkingFolder()+"WaterDepth_g.csv");

    return 0;

}
