#include "modelcreator.h"
#include "System.h"
#include "QString"
#include "fieldgenerator.h"
#include <gsl/gsl_rng.h>

ModelCreator::ModelCreator()
{

}


bool ModelCreator::Create(model_parameters mp, System *system, FieldGenerator *fieldgen)
{
    modelparameters = mp;
    TimeSeriesSet<double> SoilData;

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

    SoilData.read(path+"Soil retention params vs depth.csv");

    TimeSeriesSet<double> SoilDataCDF = SoilData.GetCummulativeDistribution();
    SoilDataCDF.write(path+"CDF.csv"); //Check CDF

    system->GetQuanTemplate(ohq_r+"main_components.json");
    system->AppendQuanTemplate(ohq_r+"unsaturated_soil_revised_model.json"); //revised version
    system->AppendQuanTemplate(ohq_r+"Well.json");
    system->AppendQuanTemplate(ohq_r+"Sewer_system.json");
    system->AppendQuanTemplate(ohq_r+"pipe_pump_tank.json");
    system->AppendQuanTemplate(ohq_r+"Pond_Plugin.json");
    system->ReadSystemSettingsTemplate(ohq_r+"settings.json");

    system->SetNumThreads(16);


    int rain_data=7; // rain data: 1: 1 yr old, 2: 1 yr new, 3: 5 yr new, 4: 3 month of 1 yr new, 5: 2 year of 1 yr new, 6: 2 days in 2022 for test

        double Simulation_start_time; // Simulation Start Date
        double Simulation_end_time; // Simulation End Date

        double Maximum_time_allowed=10*86400; // 10 Days
        double Maximum_number_of_matrix_inverstions=10*200000; // 10x

// --------------------2014-2015----------------------------------
if (rain_data==1)
    {
    Simulation_start_time=41973; // Simulation Start Date
    Simulation_end_time=42342; // Simulation End Date
    }
// --------------------2023-2024----------------------------------
else if (rain_data==2)
{
    Simulation_start_time=45230; // Simulation Start Date
    Simulation_end_time=45600; // Simulation End Date
}
// --------------------2019-2024----------------------------------
else if (rain_data==3)
{
    Simulation_start_time=43750; // Simulation Start Date
    Simulation_end_time=45600; // Simulation End Date
}
// --------------------2023-2024----------------------------------
else if (rain_data==4)
{
    Simulation_start_time=45230; // Simulation Start Date
    Simulation_end_time=45320; // Simulation End Date
}
// --------------------2022-2024----------------------------------
else if (rain_data==5)
{
    Simulation_start_time=44864; // Simulation Start Date
    Simulation_end_time=45595; // Simulation End Date
}
// --------------------2 days in 2022 for test----------------------------------
else if (rain_data==6)
{
    Simulation_start_time=44864; // Simulation Start Date
    Simulation_end_time=44866; // Simulation End Date
}
// --------------------2019 (1 days to 5 years [First run days])----------------------------------
else if (rain_data==7)
{
    Simulation_start_time=43750; // Simulation Start Date
    Simulation_end_time=43930; // Simulation End Date
}


    double dr;
    double dz;

    // Soil Blocks around gravel part

    dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    dz = mp.DepthofWell_g/mp.nz_g;

    cout<<"Soil Blocks around gravel part"<<endl;
    for (int j=0; j<mp.nz_g; j++)
    {   double actual_depth = (j+0.5)*dz+mp.DepthofWell_c;
        //calculate Ksat, n, etc.
        double Ksat = 0;
        double alpha = 0;
        double n = 0;
        double theta_s = 0;
        double theta_r = 0;

        if (!fieldgen)
        {   Ksat = SoilData["Ksat"].interpol(actual_depth,SoilDataCDF["Ksat"],mp.correlation_length_scale,false);
            alpha = SoilData["alpha"].interpol(actual_depth,SoilDataCDF["alpha"],mp.correlation_length_scale,false);
            n = SoilData["n"].interpol(actual_depth,SoilDataCDF["n"],mp.correlation_length_scale,false);
            theta_s = SoilData["theta_s"].interpol(actual_depth,SoilDataCDF["theta_s"],mp.correlation_length_scale,false);
            theta_r = SoilData["theta_r"].interpol(actual_depth,SoilDataCDF["theta_r"],mp.correlation_length_scale,false);
        }
        else
        {
            Ksat = fieldgen->interpolate("Ksat", actual_depth);
            alpha = fieldgen->interpolate("alpha", actual_depth);
            n = fieldgen->interpolate("n", actual_depth);
            theta_s = fieldgen->interpolate("theta_s", actual_depth);
            theta_r = fieldgen->interpolate("theta_r", actual_depth);
        }

        for (int i=0; i<mp.nr_g; i++)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_g + i*dr;
            double r2 = mp.rw_g + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);


            B.SetName(("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");

            if (Mode == _realization_mode::stochastic)
            {
                B.SetVal("K_sat_original",Ksat);
                B.SetVal("alpha",alpha);
                B.SetVal("n",n);
                B.SetVal("theta_sat",theta_s);
                B.SetVal("theta_res",theta_r);
            }

            else
            {
                B.SetVal("K_sat_original",SoilData["Ksat"].interpol(actual_depth));
                B.SetVal("alpha",SoilData["alpha"].interpol(actual_depth));
                B.SetVal("n",SoilData["n"].interpol(actual_depth));
                B.SetVal("theta_sat",SoilData["theta_s"].interpol(actual_depth));
                B.SetVal("theta_res",SoilData["theta_r"].interpol(actual_depth));
            }

            //B.SetVal("K_sat_original",mp.K_sat);
            //B.SetVal("K_sat_scale_factor",mp.K_o);
            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_c);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("L",mp.L);
            B.SetVal("x",-(i*dr+mp.rw_g)*2000);
            B.SetVal("y",j*dz*3000+mp.DepthofWell_c*2800);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_g);
            B.SetVal("act_Y",-actual_depth);
            system->AddBlock(B,false);
        }
    }

    // Soil Blocks under well

    dr = (mp.RadiousOfInfluence-mp.rw_uw)/mp.nr_uw;
    dz = (mp.DepthtoGroundWater-mp.DepthofWell_t)/mp.nz_uw;

    cout<<"Soil Blocks under well"<<endl;
    for (int j=0; j<mp.nz_uw; j++)
    {   double actual_depth = (j+0.5)*dz+mp.DepthofWell_c;
        //calculate Ksat, n, etc.
        double Ksat = 0;
        double alpha = 0;
        double n = 0;
        double theta_s = 0;
        double theta_r = 0;

        if (!fieldgen)
        {   Ksat = SoilData["Ksat"].interpol(actual_depth,SoilDataCDF["Ksat"],mp.correlation_length_scale,false);
            alpha = SoilData["alpha"].interpol(actual_depth,SoilDataCDF["alpha"],mp.correlation_length_scale,false);
            n = SoilData["n"].interpol(actual_depth,SoilDataCDF["n"],mp.correlation_length_scale,false);
            theta_s = SoilData["theta_s"].interpol(actual_depth,SoilDataCDF["theta_s"],mp.correlation_length_scale,false);
            theta_r = SoilData["theta_r"].interpol(actual_depth,SoilDataCDF["theta_r"],mp.correlation_length_scale,false);
        }
        else
        {
            Ksat = fieldgen->interpolate("Ksat", actual_depth);
            alpha = fieldgen->interpolate("alpha", actual_depth);
            n = fieldgen->interpolate("n", actual_depth);
            theta_s = fieldgen->interpolate("theta_s", actual_depth);
            theta_r = fieldgen->interpolate("theta_r", actual_depth);
        }

        for (int i=0; i<mp.nr_uw; i++)
        if (j*dz<mp.DepthtoGroundWater)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");
            double r1 = mp.rw_uw + i*dr;
            double r2 = mp.rw_uw + (i+1)*dr;
            double area = pi*(r2*r2-r1*r1);
            double actual_depth = (j+0.5)*dz+mp.DepthofWell_t;

            B.SetName(("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");

            if (Mode == _realization_mode::stochastic)
            {
                B.SetVal("K_sat_original",Ksat);
                B.SetVal("alpha",alpha);
                B.SetVal("n",n);
                B.SetVal("theta_sat",theta_s);
                B.SetVal("theta_res",theta_r);
            }

            else
            {
                B.SetVal("K_sat_original",SoilData["Ksat"].interpol(actual_depth));
                B.SetVal("alpha",SoilData["alpha"].interpol(actual_depth));
                B.SetVal("n",SoilData["n"].interpol(actual_depth));
                B.SetVal("theta_sat",SoilData["theta_s"].interpol(actual_depth));
                B.SetVal("theta_res",SoilData["theta_r"].interpol(actual_depth));
            }

            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("L",mp.L);
            B.SetVal("x",-(i*dr+mp.rw_uw)*2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",(i+0.5)*dr+mp.rw_uw);
            B.SetVal("act_Y",-actual_depth);
            system->AddBlock(B,false);
    }

    cout<<"Soil Blocks directly under well"<<endl;

        if (j*dz<mp.DepthtoGroundWater)
        {
            Block B;
            B.SetQuantities(system->GetMetaModel(), "Soil");

            double area = pi*pow(mp.rw_uw,2);
            double actual_depth = (j+0.5)*dz+mp.DepthofWell_t;

            B.SetName(("Soil-uw (" + QString::number(0) + "$" + QString::number(j) + ")").toStdString());
            B.SetType("Soil");

            if (Mode == _realization_mode::stochastic)
            {
                B.SetVal("K_sat_original",Ksat);
                B.SetVal("alpha",alpha);
                B.SetVal("n",n);
                B.SetVal("theta_sat",theta_s);
                B.SetVal("theta_res",theta_r);
            }

            else
            {
                B.SetVal("K_sat_original",SoilData["Ksat"].interpol(actual_depth));
                B.SetVal("alpha",SoilData["alpha"].interpol(actual_depth));
                B.SetVal("n",SoilData["n"].interpol(actual_depth));
                B.SetVal("theta_sat",SoilData["theta_s"].interpol(actual_depth));
                B.SetVal("theta_res",SoilData["theta_r"].interpol(actual_depth));
            }

            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);
            B.SetVal("theta",mp.initial_theta);
            B.SetVal("L",mp.L);
            B.SetVal("x",-mp.rw_uw*1000+2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",0);
            B.SetVal("act_Y",-actual_depth);
            system->AddBlock(B,false);
        }
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

        int j_g=mp.nz_g-1;
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

    // Well to gravel links
    cout<<"Horizontal links for well to gravel part soils"<<endl;

    dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    dz = mp.DepthofWell_g/mp.nz_g;

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

    // Graoundwater fixed head
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
    gw.SetVal("y",37000+(mp.DepthtoGroundWater)*2000);
    system->AddBlock(gw,false);

    // Groundwater links
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

    if (rain_data==1)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (1 yr).csv");
    }
    else if (rain_data==2)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (1 yr new).csv");
    }
    else if (rain_data==3)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (5 yr new).csv");
    }
    else if (rain_data==4)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (1 yr new).csv");
    }
    else if (rain_data==5)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (5 yr new).csv");
    }
    else if (rain_data==6)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (5 yr new).csv");
    }
    else if (rain_data==7)
    {
        rain.SetProperty("timeseries",path+"LA_Precipitaion (5 yr new).csv");
    }

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

    cout<<"Tracer"<<endl;
    Constituent mean_age_tracer;
    mean_age_tracer.SetQuantities(system->GetMetaModel(), "Constituent");
    mean_age_tracer.SetName("meanagetracer");
    mean_age_tracer.SetType("Constituent");
    system->AddConstituent(mean_age_tracer,false);

    cout<<"Reaction"<<endl;
    Reaction aging;
    aging.SetQuantities(system->GetMetaModel(), "Reaction");
    aging.SetName("aging");
    aging.SetType("Reaction");
    aging.SetProperty("meanagetracer:stoichiometric_constant","(1)");
    aging.SetProperty("rate_expression","(1)");
    system->AddReaction(aging,false);

    // Solve properties
    system->SetSettingsParameter("simulation_start_time",Simulation_start_time);
    system->SetSettingsParameter("simulation_end_time",Simulation_end_time);

    system->SetSettingsParameter("maximum_time_allowed",Maximum_time_allowed);
    system->SetSettingsParameter("maximum_number_of_matrix_inverstions",Maximum_number_of_matrix_inverstions);


    system->SetSystemSettings();

    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}
