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
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil_revised_model.json"); //revised version
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Sewer_system.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/pipe_pump_tank.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Pond_Plugin.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");

    double Simulation_start_time=41973; // Simulation Start Date
    double Simulation_end_time=42342; // Simulation End Date

    double dr;
    double dz;

    // Soil Blocks around gravel part

    dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    dz = mp.DepthofWell_g/mp.nz_g;

    cout<<"Soil Blocks around gravel part"<<endl;
    for (int i=0; i<mp.nr_g; i++)
        for (int j=0; j<mp.nz_g; j++)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_g + i*dr;
            double r2 = mp.rw_g + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);

            B.SetName(("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_c);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("L",mp.L);
            B.SetVal("x",-(i*dr+mp.rw_g)*2000);
            B.SetVal("y",j*dz*3000+mp.DepthofWell_c*2800);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_g);
            B.SetVal("act_Y",-(j+0.5)*dz-mp.DepthofWell_c);
            system->AddBlock(B,false);
        }

    // Soil Blocks under well

    dr = (mp.RadiousOfInfluence-mp.rw_uw)/mp.nr_uw;
    dz = (mp.DepthtoGroundWater-mp.DepthofWell_t)/mp.nz_uw;

    cout<<"Soil Blocks under well"<<endl;
    for (int i=0; i<mp.nr_uw; i++)
        for (int j=0; j<mp.nz_uw; j++)
        if (j*dz<mp.DepthtoGroundWater)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_uw + i*dr;
            double r2 = mp.rw_uw + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);

            B.SetName(("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("L",mp.L);
            B.SetVal("x",-(i*dr+mp.rw_uw)*2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_uw);
            B.SetVal("act_Y",-(j+0.5)*dz-mp.DepthofWell_t);
            system->AddBlock(B,false);
    }

    cout<<"Soil Blocks directly under well"<<endl;

    dr = (mp.RadiousOfInfluence-mp.rw_uw)/mp.nr_uw;
    dz = (mp.DepthtoGroundWater-mp.DepthofWell_t)/mp.nz_uw;

    for (int j=0; j<mp.nz_uw; j++)
        if (j*dz<mp.DepthtoGroundWater)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");

            double area = pi*pow(mp.rw_uw,2);

            B.SetName(("Soil-uw (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("L",mp.L);
            B.SetVal("x",-mp.rw_uw*1000+2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",0);
            B.SetVal("act_Y",-(j+0.5)*dz-mp.DepthofWell_t);
            system->AddBlock(B,false);
        }

    cout<<"Horizontal links for soils of gravel part"<<endl;
    for (int i=0; i<mp.nr_g; i++)
        for (int j=0; j<mp.nz_g; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");

            L.SetName(("HL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ") - Soil-g (" + QString::number(i+2) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_g));
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil-g (" + QString::number(i+2) + "$" + QString::number(j) + ")").toStdString(), false);
        }

    cout<<"Vertical links for soils of gravel part"<<endl;
    for (int i=0; i<mp.nr_g; i++)
        for (int j=0; j<mp.nz_g-1; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

            L.SetName(("VL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ") - Soil-g (" + QString::number(i+1) + "$" + QString::number(j+1)+ ")").toStdString());
            L.SetType("soil_to_soil_link");

            system->AddLink(L, ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j+1) + ")").toStdString(), false);
        }

    cout<<"Horizontal links for soils under well"<<endl;
    for (int i=0; i<mp.nr_uw; i++)
        for (int j=0; j<mp.nz_uw; j++)
            if (j*dz<mp.DepthtoGroundWater)
            {
                Link L;
                L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");

                L.SetName(("HL-Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ") - Soil-uw (" + QString::number(i+1) + "$" + QString::number(j)+ ")").toStdString());
                L.SetType("soil_to_soil_H_link");

                L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_uw));
                L.SetVal("length",dr);

                system->AddLink(L, ("Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ")").toStdString(), ("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), false);
            }

    cout<<"Vertical links for soils under well (first one)"<<endl;
    for (int i=0; i<mp.nr_uw; i++)
    {
        Link L;
        L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

        int j_g=9;
        int j_uw=0;
        L.SetName(("VL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j_g) + ") - Soil-uw (" + QString::number(i+1) + "$" + QString::number(j_uw)+ ")").toStdString());
        L.SetType("soil_to_soil_link");

        system->AddLink(L, ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j_g) + ")").toStdString(), ("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j_uw) + ")").toStdString(), false);
    }

    cout<<"Vertical links for soils under well"<<endl;
    for (int i=0; i<mp.nr_uw+1; i++)
        for (int j=0; j<mp.nz_uw; j++)
            if (j*dz<mp.DepthtoGroundWater)
            {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

            L.SetName(("VL-Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ") - Soil-uw (" + QString::number(i) + "$" + QString::number(j+1)+ ")").toStdString());
            L.SetType("soil_to_soil_link");

            system->AddLink(L, ("Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ")").toStdString(), ("Soil-uw (" + QString::number(i) + "$" + QString::number(j+1) + ")").toStdString(), false);
    }


    cout<<"Well_c"<<endl;
    Block well_c;
    well_c.SetQuantities(system->GetMetaModel(), "Well_aggregate");
    well_c.SetName("Well_c");
    well_c.SetType("Well_aggregate");
    well_c.SetVal("_height",mp.DepthofWell_c*2000);
    well_c.SetVal("_width",mp.rw_c*1000);
    well_c.SetVal("bottom_elevation",-mp.DepthofWell_c);
    well_c.SetVal("diameter",mp.rw_c*2);
    well_c.SetVal("depth",mp.depth_w_c);
    well_c.SetVal("porosity",mp.porosity_c);
    well_c.SetVal("x",-mp.rw_c*1000+2000);
    well_c.SetVal("y",mp.DepthofWell_c*200);
    system->AddBlock(well_c,false);

    cout<<"Well_g"<<endl;
    Block well_g;
    well_g.SetQuantities(system->GetMetaModel(), "Well_aggregate");
    well_g.SetName("Well_g");
    well_g.SetType("Well_aggregate");
    well_g.SetVal("_height",mp.DepthofWell_g*3200);
    well_g.SetVal("_width",mp.rw_g*1000);
    well_g.SetVal("bottom_elevation",-(mp.DepthofWell_c+mp.DepthofWell_g));
    well_g.SetVal("diameter",mp.rw_g*2);
    well_g.SetVal("depth",mp.depth_w_g);
    well_g.SetVal("porosity",mp.porosity_g);
    well_g.SetVal("x",-mp.rw_g*1000+2000);
    well_g.SetVal("y",(mp.DepthofWell_c+mp.DepthofWell_g)*1000);
    system->AddBlock(well_g,false);

    cout<<"Well to well overflow"<<endl;
    Link well_to_well_overflow;
    well_to_well_overflow.SetQuantities(system->GetMetaModel(), "Sewer_pipe");
    well_to_well_overflow.SetName("Well_to_well_overflow");
    well_to_well_overflow.SetType("Sewer_pipe");
    well_to_well_overflow.SetVal("ManningCoeff",mp.ManningCoeff_of);
    well_to_well_overflow.SetVal("diameter",mp.diameter_of);
    well_to_well_overflow.SetVal("length",mp.length_of);
    well_to_well_overflow.SetVal("start_elevation",mp.start_elevation_of);
    well_to_well_overflow.SetVal("end_elevation",mp.end_elevation_of);
    system->AddLink(well_to_well_overflow,"Well_c","Well_g",false);

    cout<<"Junction"<<endl;
    Block junction;
    junction.SetQuantities(system->GetMetaModel(), "junction_elastic");
    junction.SetName("Junction_elastic");
    junction.SetType("junction_elastic");
    junction.SetVal("_height",1000);
    junction.SetVal("_width",1000);
    junction.SetVal("x",3000);
    junction.SetVal("y",mp.DepthofWell_c*2000+1000);
    junction.SetVal("elevation",mp.elevation_j); // Should be calculated
    system->AddBlock(junction,false);

    cout<<"Well to junction"<<endl;
    Link well_to_junction;
    well_to_junction.SetQuantities(system->GetMetaModel(), "darcy_connector");
    well_to_junction.SetName("Well_to_junction");
    well_to_junction.SetType("darcy_connector");
    system->AddLink(well_to_junction,"Well_c","Junction_elastic",false);

    cout<<"Junction to well"<<endl;
    Link junction_to_well;
    junction_to_well.SetQuantities(system->GetMetaModel(), "darcy_connector");
    junction_to_well.SetName("Junction_to_well");
    junction_to_well.SetType("darcy_connector");
    system->AddLink(junction_to_well,"Junction_elastic","Well_g",false);

    cout<<"Horizontal links for well to gravel part soils"<<endl;
    for (int j=0; j<mp.nz_g; j++)
        if (j*dz<mp.DepthofWell_t)
        {
                Link L;
                L.SetQuantities(system->GetMetaModel(), "Well2soil horizontal link");

                int i_g=1;
                L.SetName(("HL_Well_g - Soil-g (" + QString::number(i_g) + "$" + QString::number(j)+ ")").toStdString());
                L.SetType("Well2soil horizontal link");

                L.SetVal("length",dr/2);

                system->AddLink(L, "Well_g", ("Soil-g (" + QString::number(i_g) + "$" + QString::number(j) + ")").toStdString(), false);
        }

    cout<<"Vertical links for well to under well soils"<<endl;

        Link L1;
        L1.SetQuantities(system->GetMetaModel(), "Well2soil vertical link");

        int i_uw=0;
        int j_uw=0;
        L1.SetName(("VL_Well_g - Soil-uw (" + QString::number(i_uw) + "$" + QString::number(j_uw)+ ")").toStdString());
        L1.SetType("Well2soil vertical link");

        system->AddLink(L1, "Well_g", ("Soil-uw (" + QString::number(i_uw) + "$" + QString::number(j_uw) + ")").toStdString(), false);

    cout<<"Groundwater"<<endl;
    Block gw;
    gw.SetQuantities(system->GetMetaModel(), "fixed_head");
    gw.SetName("Ground Water");
    gw.SetType("fixed_head");
    gw.SetVal("_height",500);
    gw.SetVal("_width",mp.RadiousOfInfluence*1000);
    gw.SetVal("head",-mp.DepthtoGroundWater);
    gw.SetVal("Storage",100000);
    gw.SetVal("x",-mp.nr_uw*1000);
    gw.SetVal("y",mp.DepthtoGroundWater*2200);
    system->AddBlock(gw,false);

    cout<<"Soil to Groundwater"<<endl;
    for (int i=0; i<mp.nr_uw+1; i++)
    {
        Link L2;
        L2.SetQuantities(system->GetMetaModel(), "soil_to_fixedhead_link");
        L2.SetName(("Soil to Groundwater (" + QString::number(i) + ")").toStdString());
        L2.SetType("soil_to_fixedhead_link");

        system->AddLink(L2, ("Soil-uw (" + QString::number(i) + "$" + QString::number(mp.nz_uw-1) + ")").toStdString(), "Ground Water", false);
    }

    cout<<"Rain"<<endl;
    Source rain;
    rain.SetQuantities(system->GetMetaModel(), "Precipitation");
    rain.SetType("Precipitation");
    rain.SetName("Rain");
    rain.SetVal("_height",3000);
    rain.SetVal("_width",3000);
    rain.SetVal("x",-5000);
    rain.SetVal("y",-500);
    rain.SetProperty("timeseries","/home/behzad/Projects/VN Drywell_Models/LA_Precipitaion (1 yr).csv"); // Should be modified according to the users directory
    system->AddSource(rain, false);

    cout<<"Catchment"<<endl;
    Block catchment;
    catchment.SetQuantities(system->GetMetaModel(), "Catchment");
    catchment.SetName("Catchment");
    catchment.SetType("Catchment");
    catchment.SetVal("_height",3000);
    catchment.SetVal("_width",3000);
    catchment.SetVal("x",-5000);
    catchment.SetVal("y",0);
    catchment.SetVal("ManningCoeff",mp.ManningCoeff_cm);
    catchment.SetVal("Slope",mp.Slope_cm);
    catchment.SetVal("area",mp.area_cm);
    catchment.SetVal("Width",mp.Width_cm);
    catchment.SetVal("y",0);
    system->AddBlock(catchment,false);
    system->block("Catchment")->SetProperty("Precipitation","Rain");
    cout<<"Catchment to well link"<<endl;

    Link L3;
    L3.SetQuantities(system->GetMetaModel(), "Surface water to well");
    L3.SetName("Catchment to well");
    L3.SetType("Surface water to well");
    L3.SetVal("ManningCoeff",mp.ManningCoeff_cmw);
    L3.SetVal("length",mp.length_cmw);
    system->AddLink(L3, "Catchment", "Well_c", false);


    system->SetSettingsParameter("simulation_start_time",Simulation_start_time);
    system->SetSettingsParameter("simulation_end_time",Simulation_end_time);

    system->SetSystemSettings();

    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}




/*
    double dr = (mp.RadiousOfInfluence-mp.rw_c_t)/mp.nr_c;
    double dz = mp.DepthofWell_c/mp.nz_c;

    // Soil Blocks around concrete part
    cout<<"Soil Blocks around concrete part"<<endl;
    for (int i=0; i<mp.nr_c; i++)
        for (int j=0; j<mp.nz_c; j++)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_c_t + i*dr;
            double r2 = mp.rw_c_t + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);

            B.SetName(("Soil-c (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");
            B.SetVal("K_sat_original",mp.K_sat);
            B.SetVal("alpha",mp.alpha);
            B.SetVal("area",area);
            B.SetVal("_width",dr*300);
            B.SetVal("_height",dr*300);
            B.SetVal("bottom_elevation",-(j+1)*dz);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("theta_sat",mp.theta_sat);
            B.SetVal("theta_res",mp.theta_r);
            B.SetVal("x",-(i*dr+mp.rw_c_t)*2000);
            B.SetVal("y",1000+(j*dz)*2000);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_c_t);
            B.SetVal("act_Y",-(j+0.5)*dz);
            system->AddBlock(B,false);
        }
    */

/*
    cout<<"Horizontal links for soils of concrete part"<<endl;
    for (int i=0; mp.nr_c<i<mp.RadiousOfInfluence; i++)
        for (int j=0; mp.nz_c<j<mp.DepthtoGroundWater; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");


            L.SetName(("Soil-c (" + QString::number(i+1) + "$" + QString::number(j) + ") - Soil-c (" + QString::number(i+2) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_g));
            L.SetVal("length",dr);

            system->AddLink(L, ("Soil-c (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(), ("Soil-c (" + QString::number(i+2) + "$" + QString::number(j) + ")").toStdString(), false);
        }
    */
