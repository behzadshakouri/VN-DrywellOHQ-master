#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"


int main(int argc, char *argv[])
{

    model_parameters mp;

    System *system=new System();
    ModelCreator ModCreate;
    cout<<"Creating model ..." <<endl;
    ModCreate.Create(mp,system);
    cout<<"Creating model done..." <<endl;

    system->SetWorkingFolder("/mnt/3rd900/Projects/VN Drywell_Models/"); // Should be modified according to the users directory
    system->SetSilent(false);
    cout<<"Saving"<<endl;
    system->SavetoScriptFile("/mnt/3rd900/Projects/VN Drywell_Models/CreatedModel.ohq"); // Should be modified according to the users directory

    cout<<"Solving ..."<<endl;
    system->Solve();
    system->SavetoJson("/mnt/3rd900/Projects/VN Drywell_Models/Model.json",system->addedtemplates, true, true );
    cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() +"'"<<endl;
    CTimeSeriesSet<double> uniformoutput_LR = system->GetOutputs().make_uniform(1);
    CTimeSeriesSet<double> uniformoutput_HR = system->GetOutputs().make_uniform(0.1);
    uniformoutput_HR.writetofile(system->GetWorkingFolder() + system->OutputFileName());
    cout<<"Getting results into grid"<<endl;
    ResultGrid resgrid(uniformoutput_LR,"theta",system);
    cout<<"Writing VTPs"<<endl;
    resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"moisture.vtp");

    vector<string> well_block; well_block.push_back("Well");

    ResultGrid well_depth = ResultGrid(uniformoutput_HR,well_block,"depth" );
    well_depth.Sum().writefile(system->GetWorkingFolder()+"WaterDepth.csv");

    return 0;

}
