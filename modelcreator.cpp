#include "modelcreator.h"
#include "System.h"
#include "QString"

#include <gsl/gsl_rng.h>

ModelCreator::ModelCreator()
{

}


bool ModelCreator::Create(model_parameters mp, System *system)
{
    system->GetQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/main_components.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");

    double dr = (mp.RadiousOfInfluence-mp.rw)/mp.nr;
    double dz = mp.DepthtoGroundWater/mp.nz;

    cout<<"Main Soil Blocks"<<endl;
    for (int i=0; i<mp.nr; i++)
        for (int j=0; j<mp.nz; j++)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw + i*dr;
            double r2 = mp.rw + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);

            B.SetName(("Soil (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*3000);
            B.SetVal("_height",dr*3000);
            B.SetVal("bottom_elevation",-(j+1)*dz);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("x",(i*dr+mp.rw)*4000);
            B.SetVal("y",(j*dz)*4000);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw);
            B.SetVal("act_Y",-(j+0.5)*dz);
            system->AddBlock(B,false);
        }

    cout<<"Soil Blocks under well"<<endl;
    for (int j=0; j<mp.nz; j++)
    {
        if (j*dz>mp.DepthofWell)
        {   Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");

            double area = pi*(mp.rw*mp.rw);

            B.SetName(("Soil (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*3000);
            B.SetVal("_height",dr*3000);
            B.SetVal("bottom_elevation",-(j+1)*dz);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("x",-mp.rw*4000);
            B.SetVal("y",(j*dz)*4000);
            B.SetVal("act_X",0);
            B.SetVal("act_Y",-(j+0.5)*dz);
            system->AddBlock(B,false);
        }
    }


    cout<<"Horizontal links"<<endl;
    for (int i=0; i<mp.nr-1; i++)
        for (int j=0; j<mp.nz; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");


            L.SetName(("Soil (" + QString::number(i+1) + "$" + QString::number(j) + "-" + QString::number(i+2) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw));
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(i+2) + "$" + QString::number(j) + ")").toStdString(), false);
        }

    cout<<"Horizontal links for central soils"<<endl;
    for (int j=0; j<mp.nz; j++)
    {

        if (j*dz>mp.DepthofWell)
        {   Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");


            L.SetName(("Soil (" + QString::number(0) + "$" + QString::number(j) + "-" + QString::number(1) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*mp.rw);
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(1) + "$" + QString::number(j) + ")").toStdString(), false);
        }
    }
    int well_layer=0;
    cout<<"Vertical links for central soils"<<endl;
    for (int j=0; j<mp.nz-1; j++)
    {

        if (j*dz>mp.DepthofWell)
        {   Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");
            if (well_layer==0) well_layer = j;

            L.SetName(("Soil (" + QString::number(0) + "$" + QString::number(j) + "-" + QString::number(0) + "$" + QString::number(j+1)+ ")").toStdString());
            L.SetType("soil_to_soil_link");

            system->AddLink(L, ("Soil (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(0) + "$" + QString::number(j+1) + ")").toStdString(), false);
        }
    }
    cout<<"Vertical links"<<endl;
    for (int i=0; i<mp.nr; i++)
        for (int j=0; j<mp.nz-1; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");


            L.SetName(("Soil (" + QString::number(i+1) + "$" + QString::number(j) + "-" + QString::number(i) + "$" + QString::number(j+1)+ ")").toStdString());
            L.SetType("soil_to_soil_link");

            //L.SetVal("length",dr);

            system->AddLink(L, ("Soil (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(i+1) + "$" + QString::number(j+1) + ")").toStdString(), false);
        }
    cout<<"Well"<<endl;
    Block well;
    well.SetQuantities(system->GetMetaModel(), "Well_aggregate");
    well.SetName("Well");
    well.SetType("Well_aggregate");
    well.SetVal("_height",mp.DepthofWell*4000);
    well.SetVal("_width",mp.rw*4000);
    well.SetVal("bottom_elevation",-mp.DepthofWell);
    well.SetVal("diameter",mp.rw*2);
    well.SetVal("depth",0.5);
    well.SetVal("porosity",1);
    well.SetVal("x",-mp.rw*4000);
    well.SetVal("y",0);
    system->AddBlock(well,false);


    cout<<"Well to soil"<<endl;

    for (int j=0; j<mp.nz; j++)
    {
        if (j*dz<mp.DepthofWell)
        {   Link L;
            L.SetQuantities(system->GetMetaModel(), "Well2soil horizontal link");
            L.SetName(("Well to Soil (" + QString::number(j) + ")").toStdString());
            L.SetType("Well2soil horizontal link");
            L.SetVal("length",dr/2);

            system->AddLink(L, "Well", ("Soil (" + QString::number(1) + "$" + QString::number(j) + ")").toStdString(), false);
        }
    }

    cout<<"Well to bottom"<<endl;
    Link well_to_bottom;
    well_to_bottom.SetQuantities(system->GetMetaModel(), "Well2soil vertical link");
    well_to_bottom.SetName("Well_to_bottom");
    well_to_bottom.SetType("Well2soil vertical link");
    system->AddLink(well_to_bottom,"Well",("Soil (" + QString::number(0) + "$" + QString::number(well_layer) + ")").toStdString(), false);

    cout<<"Groundwater"<<endl;
    Block gw;
    gw.SetQuantities(system->GetMetaModel(), "fixed_head");
    gw.SetName("Ground Water");
    gw.SetType("fixed_head");
    gw.SetVal("_height",200);
    gw.SetVal("_width",mp.RadiousOfInfluence*4000);
    gw.SetVal("head",-mp.DepthtoGroundWater);
    gw.SetVal("Storage",100000);
    gw.SetVal("x",0);
    gw.SetVal("y",mp.DepthtoGroundWater*4000+400);
    system->AddBlock(gw,false);

    cout<<"Soil to Groundwater"<<endl;
    for (int i=0; i<mp.nr+1; i++)
    {
        Link L;
        L.SetQuantities(system->GetMetaModel(), "soil_to_fixedhead_link");
        L.SetName(("Soil to Groundwater (" + QString::number(i) + ")").toStdString());
        L.SetType("soil_to_fixedhead_link");

        system->AddLink(L, ("Soil (" + QString::number(i) + "$" + QString::number(mp.nz-1) + ")").toStdString(), "Ground Water", false);

    }

    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}
