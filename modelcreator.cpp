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

    // Soil blocks around concrete part
    double dr = (mp.RadiousOfInfluence-mp.rw_c_t)/mp.nr;
    double dz = mp.DepthtoGroundWater/mp.nz;

    cout<<"Main Soil Blocks"<<endl;
    for (int i=0; i<mp.nr; i++)
        for (int j=0; j<mp.nz; j++)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_c_t + i*dr;
            double r2 = mp.rw_c_t + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);

            B.SetName(("Soil (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*200);
            B.SetVal("_height",dr*200);
            B.SetVal("bottom_elevation",-(j+1)*dz);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("x",-(i*dr+mp.rw_c_t)*1000);
            B.SetVal("y",(j*dz)*1000);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_c_t);
            B.SetVal("act_Y",-(j+0.5)*dz);
            system->AddBlock(B,false);
        }

    cout<<"Soil Blocks under wells"<<endl;
    for (int j=0; j<mp.nz; j++)
    {
        if (j*dz>mp.DepthofWell_c)
        {   Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");

            double area = pi*(mp.rw_c_t*mp.rw_c_t);

            B.SetName(("Soil (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*200);
            B.SetVal("_height",dr*200);
            B.SetVal("bottom_elevation",-(j+1)*dz);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("x",-mp.rw_c_t*1000);
            B.SetVal("y",(j*dz)*1000);
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

            L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_c_t));
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(i+2) + "$" + QString::number(j) + ")").toStdString(), false);
        }

    cout<<"Horizontal links for central soils"<<endl;
    for (int j=0; j<mp.nz; j++)
    {

        if (j*dz>mp.DepthofWell_c)
        {   Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");


            L.SetName(("Soil (" + QString::number(0) + "$" + QString::number(j) + "-" + QString::number(1) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*mp.rw_c_t);
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString(), ("Soil (" + QString::number(1) + "$" + QString::number(j) + ")").toStdString(), false);
        }
    }
    int well_layer=0;
    cout<<"Vertical links for central soils"<<endl;
    for (int j=0; j<mp.nz-1; j++)
    {

        if (j*dz>mp.DepthofWell_c)
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

    cout<<"Well_c"<<endl;
    Block well_c;
    well_c.SetQuantities(system->GetMetaModel(), "Well_aggregate");
    well_c.SetName("Well_c");
    well_c.SetType("Well_aggregate");
    well_c.SetVal("_height",mp.DepthofWell_c*1000);
    well_c.SetVal("_width",mp.rw_c*1000);
    well_c.SetVal("bottom_elevation",-mp.DepthofWell_c);
    well_c.SetVal("diameter",mp.rw_c*2);
    well_c.SetVal("depth",0.5);
    well_c.SetVal("porosity",mp.porosity_c);
    well_c.SetVal("x",-mp.rw_c*1000);
    well_c.SetVal("y",mp.DepthofWell_c*1000);
    system->AddBlock(well_c,false);

    cout<<"Well_g1"<<endl;
    Block well_g;
    well_g.SetQuantities(system->GetMetaModel(), "Well_aggregate");
    well_g.SetName("Well_g1");
    well_g.SetType("Well_aggregate");
    well_g.SetVal("_height",mp.DepthofWell_g*1000);
    well_g.SetVal("_width",mp.rw_g*1000);
    well_g.SetVal("bottom_elevation",-(mp.DepthofWell_c+mp.DepthofWell_g));
    well_g.SetVal("diameter",mp.rw_g*2);
    well_g.SetVal("depth",0.5);
    well_g.SetVal("porosity",mp.porosity_g);
    well_g.SetVal("x",-mp.rw_g*1000);
    well_g.SetVal("y",(mp.DepthofWell_c+mp.DepthofWell_g)*2*1000);
    system->AddBlock(well_g,false);

    cout<<"Well to soil"<<endl;

    for (int j=0; j<mp.nz; j++)
    {
        if (j*dz<mp.DepthofWell_t)
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
    gw.SetVal("_width",mp.RadiousOfInfluence*1000);
    gw.SetVal("head",-mp.DepthtoGroundWater);
    gw.SetVal("Storage",100000);
    gw.SetVal("x",0);
    gw.SetVal("y",mp.DepthtoGroundWater*1000+100);
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
