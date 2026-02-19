#include "modelcreator.h"
#include "System.h"
#include "QString"
#include "fieldgenerator.h"
#include <gsl/gsl_rng.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>
#include <string>

// -----------------------------
// Local helpers (no header changes elsewhere)
// -----------------------------
static inline bool file_exists_mc(const std::string& p)
{
    std::ifstream f(p);
    return (bool)f;
}

// Reads FIRST time-slice profile from a tidy CSV like:
// borehole,time_str,time_excel,depth_m,soil_moisture_pct
// Returns depth_m (meters) and theta (fraction 0..1)
static bool read_ert_tidy_first_profile_mc(
    const std::string& path,
    std::vector<double>& depth_m,
    std::vector<double>& theta_frac,
    double depth_offset_m = 0.0,
    double time_tol = 1e-9)
{
    depth_m.clear();
    theta_frac.clear();

    std::ifstream f(path);
    if (!f) return false;

    struct Row { double t; double d; double th; };
    std::vector<Row> rows;
    rows.reserve(256);

    std::string line;
    bool first = true;

    while (std::getline(f, line))
    {
        if (line.empty()) continue;

        if (first)
        {
            first = false;
            bool has_alpha = std::any_of(line.begin(), line.end(), [](unsigned char c){
                return std::isalpha(c);
            });
            if (has_alpha) continue; // skip header
        }

        // borehole,time_str,time_excel,depth_m,soil_moisture_pct
        std::stringstream ss(line);
        std::string c0, c1, c2, c3, c4;
        if (!std::getline(ss, c0, ',')) continue;
        if (!std::getline(ss, c1, ',')) continue;
        if (!std::getline(ss, c2, ',')) continue;
        if (!std::getline(ss, c3, ',')) continue;
        if (!std::getline(ss, c4, ',')) continue;

        try {
            double t_excel = std::stod(c2);
            double d       = std::stod(c3) + depth_offset_m;

            // tidy has theta in percent (0..100) -> model needs 0..1
            double sm_pct  = std::stod(c4);
            double th      = sm_pct;

            // "check" behavior: if values look like percent, scale; otherwise keep
            if (th > 1.5) th /= 100.0;

            // hard clamp to [0,1]
            if (th < 0.0) th = 0.0;
            if (th > 1.0) th = 1.0;

            rows.push_back({t_excel, d, th});
        } catch (...) {
            // skip bad lines
        }
    }

    if (rows.empty()) return false;

    // Find earliest time
    double t0 = rows.front().t;
    for (const auto& r : rows) t0 = std::min(t0, r.t);

    // Keep only rows at earliest time slice
    for (const auto& r : rows)
    {
        if (std::abs(r.t - t0) <= time_tol)
        {
            depth_m.push_back(r.d);
            theta_frac.push_back(r.th);
        }
    }

    if (depth_m.empty()) return false;

    // Sort by depth
    if (depth_m.size() >= 2)
    {
        std::vector<size_t> idx(depth_m.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

        std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j){
            return depth_m[i] < depth_m[j];
        });

        std::vector<double> d2, t2;
        d2.reserve(depth_m.size());
        t2.reserve(theta_frac.size());

        for (size_t k = 0; k < idx.size(); ++k) {
            d2.push_back(depth_m[idx[k]]);
            t2.push_back(theta_frac[idx[k]]);
        }

        depth_m.swap(d2);
        theta_frac.swap(t2);
    }

    return true;
}

static double linear_interp_or_nan_mc(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double xq)
{
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    if (x.size() < 2 || x.size() != y.size()) return NaN;
    if (xq < x.front() || xq > x.back()) return NaN;

    auto it = std::lower_bound(x.begin(), x.end(), xq);
    if (it == x.begin()) return y.front();
    if (it == x.end()) return y.back();

    size_t hi = (size_t)std::distance(x.begin(), it);
    size_t lo = hi - 1;

    const double x0 = x[lo], x1 = x[hi];
    const double y0 = y[lo], y1 = y[hi];

    const double dx = x1 - x0;
    if (std::abs(dx) < 1e-30) return y0;

    const double w = (xq - x0) / dx;
    return y0 + w * (y1 - y0);
}

struct ObsBH_MC {
    std::string name;
    double r_m = 0.0;
    std::vector<double> depth;
    std::vector<double> theta;
    double depth_offset_m = 0.0;
};

// IDW in radius using available boreholes at requested depth
static double theta_idw_r_from_bhs_mc(
    const std::vector<ObsBH_MC>& bhs,
    double r_m,
    double depth_m,
    double fallback_theta)
{
    const double eps = 0.05;  // m
    const double p   = 2.0;

    double num = 0.0, den = 0.0;

    for (const auto& bh : bhs)
    {
        double th = std::numeric_limits<double>::quiet_NaN();

        if (bh.depth.size() >= 2) th = linear_interp_or_nan_mc(bh.depth, bh.theta, depth_m);
        else if (bh.depth.size() == 1 && std::abs(bh.depth[0] - depth_m) < 1e-6) th = bh.theta[0];

        if (!std::isfinite(th)) continue;

        double w = 1.0 / std::pow(std::abs(r_m - bh.r_m) + eps, p);
        num += w * th;
        den += w;
    }

    if (den > 0.0) return num / den;
    return fallback_theta;
}

// Average across available boreholes (no r dependence)
static double theta_r_avg_profile_mc(
    const std::vector<ObsBH_MC>& bhs,
    double depth_m,
    double fallback_theta)
{
    double sum = 0.0;
    int    cnt = 0;

    for (const auto& bh : bhs)
    {
        double th = std::numeric_limits<double>::quiet_NaN();

        if (bh.depth.size() >= 2) th = linear_interp_or_nan_mc(bh.depth, bh.theta, depth_m);
        else if (bh.depth.size() == 1 && std::abs(bh.depth[0] - depth_m) < 1e-6) th = bh.theta[0];

        if (!std::isfinite(th)) continue;

        sum += th;
        cnt++;
    }

    if (cnt > 0) return sum / (double)cnt;
    return fallback_theta;
}

// -----------------------------
// CLI parsing (declared in modelcreator.h)
// -----------------------------
static bool starts_with_mc(const std::string& s, const std::string& pfx)
{
    return s.size() >= pfx.size() && s.compare(0, pfx.size(), pfx) == 0;
}

// Accepts: --ksat-scale 0.5   or   --ksat-scale=0.5
double parse_ksat_scale(int argc, char* argv[], double default_val)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i] ? std::string(argv[i]) : std::string();

        std::string v;
        const std::string key1 = "--ksat-scale=";
        const std::string key2 = "--ksat-scale";

        if (starts_with_mc(a, key1))
        {
            v = a.substr(key1.size());
        }
        else if (a == key2)
        {
            if (i + 1 < argc) v = argv[i + 1] ? std::string(argv[i + 1]) : std::string();
        }
        else
        {
            continue;
        }

        try {
            double x = std::stod(v);
            if (std::isfinite(x) && x > 0.0) return x;
        } catch (...) {
            return default_val;
        }
    }

    return default_val;
}

// -----------------------------
// ModelCreator
// -----------------------------
ModelCreator::ModelCreator()
{
}

bool ModelCreator::Create(model_parameters mp,
                          System* system,
                          FieldGenerator* fieldgen,
                          const RainConfig& raincfg,
                          const SimulationConfig& simcfg)
{
    modelparameters = mp;
    TimeSeriesSet<double> SoilData;

#ifdef Behzad
    const std::string path="/home/behzad/Projects/VN Drywell_Models/";
    const std::string ohq_r="/home/behzad/Projects/OpenHydroQual/resources/";
#elif PowerEdge
    const std::string path="/mnt/3rd900/Projects/VN Drywell_Models/";
    const std::string ohq_r="/mnt/3rd900/Projects/OpenHydroQual/resources/";
#elif Arash
    const std::string path="/home/arash/Projects/VN Drywell_Models/";
    const std::string ohq_r="/home/arash/Projects/OpenHydroQual/resources/";
#elif SligoCreek
    const std::string path="/media/arash/E/Projects/VN Drywell_Models/";
    const std::string ohq_r="/media/arash/E/Projects/OpenHydroQual/resources/";
#endif

    SoilData.read(path+"Soil retention params vs depth.csv");

    TimeSeriesSet<double> SoilDataCDF = SoilData.GetCummulativeDistribution();
    SoilDataCDF.write(path+"CDF.csv"); // Check CDF

    system->GetQuanTemplate(ohq_r+"main_components.json");
    system->AppendQuanTemplate(ohq_r+"unsaturated_soil_revised_model.json"); // revised version
    system->AppendQuanTemplate(ohq_r+"Well.json");
    system->AppendQuanTemplate(ohq_r+"Sewer_system.json");
    system->AppendQuanTemplate(ohq_r+"pipe_pump_tank.json");
    system->AppendQuanTemplate(ohq_r+"Pond_Plugin.json");
    system->ReadSystemSettingsTemplate(ohq_r+"settings.json");

    system->SetNumThreads(16);

    // ==================================================================
    // Ksat scale factor (controlled from main via simcfg.KsatScaleFactor)
    // unsaturated_soil_revised_model.json should use:
    //   K_sat = K_sat_original * K_sat_scale_factor
    // ==================================================================
    double KsatScale = simcfg.KsatScaleFactor;
    if (!std::isfinite(KsatScale) || KsatScale <= 0.0) KsatScale = 1.0;

    // ==================================================================
    // OPTIONAL: Assign initial theta from ERT tidy observed profiles
    // Uses FIRST time-slice in each tidy file and maps pct->fraction
    // Honors:
    //   ERTInitialThetaMode  (0=IDW through r, 1=avg through r)
    //   ERTInitialThetaWhich (0=both, 3=ERT-3 only, 5=ERT-5 only)
    // ==================================================================
    std::vector<ObsBH_MC> initBHs;

    if (UseERTInitialTheta)
    {
        auto addBH = [&](const std::string& name, double r_m, const std::string& csv, double depth_offset_m)
        {
            ObsBH_MC bh;
            bh.name = name;
            bh.r_m  = r_m;
            bh.depth_offset_m = depth_offset_m;

            std::vector<double> d, t;
            if (read_ert_tidy_first_profile_mc(csv, d, t, depth_offset_m))
            {
                bh.depth = std::move(d);
                bh.theta = std::move(t);
                initBHs.push_back(std::move(bh));
                std::cout << "ERT-IC: loaded FIRST time-slice " << name
                          << " (" << csv << ") n=" << initBHs.back().depth.size() << "\n";
            }
            else
            {
                std::cout << "ERT-IC: WARNING could not read " << csv << " for " << name << "\n";
            }
        };

        const bool want3 = (ERTInitialThetaWhich == 0 || ERTInitialThetaWhich == 3);
        const bool want5 = (ERTInitialThetaWhich == 0 || ERTInitialThetaWhich == 5);

        if (want3) addBH("ERT-3", 6.7979, system->GetWorkingFolder() + "ERT-3_tidy.csv", 0.0);
        if (want5) addBH("ERT-5", 4.3974, system->GetWorkingFolder() + "ERT-5_tidy.csv", 0.0);

        if (initBHs.empty())
            std::cout << "ERT-IC: No valid borehole tidy loaded -> will fallback to mp.initial_theta\n";
    }

    auto theta_from_ert = [&](double r_m, double depth_m, double fallback)->double
    {
        if (!UseERTInitialTheta || initBHs.empty()) return fallback;

        if (ERTInitialThetaMode == 1)
            return theta_r_avg_profile_mc(initBHs, depth_m, fallback); // avg through r
        else
            return theta_idw_r_from_bhs_mc(initBHs, r_m, depth_m, fallback); // IDW in r
    };

    // ==================================================================
    // RAINFALL DATABASE
    // ==================================================================
    struct RainDataset {
        std::string file;
        double start;
        double end;
    };

    std::map<int, RainDataset> rainDB = {
        {1, { "LA_Precipitaion (1 yr).csv",                     41973, 42342 }},
        {2, { "LA_Precipitaion (1 yr new).csv",                 45230, 45600 }},
        {3, { "LA_Precipitaion (5 yr new).csv",                 44864, 45595 }},
        {4, { "Pacoima Spreading Grounds_rainfall data.csv",     43466, 45292 }},
        {5, { "Synthetic_rain_flow.csv",                        45763, 45764 }},
    };

    if (!rainDB.count(raincfg.rain_data)) {
        std::cerr << "Invalid rain_data=" << raincfg.rain_data << "\n";
        return false;
    }

    auto RC = rainDB[raincfg.rain_data];

    double Simulation_start_time =
        (simcfg.start_time > 0 ? simcfg.start_time : RC.start);

    double Simulation_end_time =
        (simcfg.end_time > 0 ? simcfg.end_time : RC.end);

    std::string rain_file = path + RC.file;

    // ==================================================================
    // Model blocks, links, ...
    // ==================================================================
    double dr;
    double dz;

    int nz_uw_n = (raincfg.rain_data == 5) ? mp.nz_uw_n : mp.nz_uw;

    // Soil Blocks around gravel part
    dr = (mp.RadiousOfInfluence-mp.rw_g)/mp.nr_g;
    dz = mp.DepthofWell_g/mp.nz_g;

    std::cout<<"Soil Blocks around gravel part"<<std::endl;
    for (int j=0; j<mp.nz_g; j++)
    {
        double actual_depth = (j+0.5)*dz+mp.DepthofWell_c;

        double Ksat = 0;
        double alpha = 0;
        double n = 0;
        double theta_s = 0;
        double theta_r = 0;

        if (!fieldgen)
        {
            Ksat = SoilData["Ksat"].interpol(actual_depth,SoilDataCDF["Ksat"],mp.correlation_length_scale,false);
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

            // NEW: Ksat scale factor (used by unsaturated_soil_revised_model.json)
            B.SetVal("K_sat_scale_factor", KsatScale);

            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_c);
            B.SetVal("depth",dz);

            // Geometry
            double actX = (i+0.5)*dr + mp.rw_g;  // radius (m)
            double actY = -actual_depth;         // depth stored negative

            B.SetVal("x",-(i*dr+mp.rw_g)*2000);
            B.SetVal("y",j*dz*3000+mp.DepthofWell_c*2800);
            B.SetVal("act_X",actX);
            B.SetVal("act_Y",actY);

            // INITIAL THETA
            double theta_init = theta_from_ert(actX, actual_depth, mp.initial_theta);
            B.SetVal("theta",theta_init);

            B.SetVal("L",mp.L);

            system->AddBlock(B,false);
        }
    }

    // Soil Blocks under well
    dr = (mp.RadiousOfInfluence-mp.rw_uw)/mp.nr_uw;
    dz = (mp.DepthtoGroundWater-mp.DepthofWell_t)/mp.nz_uw;

    std::cout<<"Soil Blocks under well"<<std::endl;
    for (int j = 0; j < nz_uw_n; j++)
    {
        double actual_depth_outer = (j+0.5)*dz+mp.DepthofWell_c;

        double Ksat = 0;
        double alpha = 0;
        double n = 0;
        double theta_s = 0;
        double theta_r = 0;

        if (!fieldgen)
        {
            Ksat = SoilData["Ksat"].interpol(actual_depth_outer,SoilDataCDF["Ksat"],mp.correlation_length_scale,false);
            alpha = SoilData["alpha"].interpol(actual_depth_outer,SoilDataCDF["alpha"],mp.correlation_length_scale,false);
            n = SoilData["n"].interpol(actual_depth_outer,SoilDataCDF["n"],mp.correlation_length_scale,false);
            theta_s = SoilData["theta_s"].interpol(actual_depth_outer,SoilDataCDF["theta_s"],mp.correlation_length_scale,false);
            theta_r = SoilData["theta_r"].interpol(actual_depth_outer,SoilDataCDF["theta_r"],mp.correlation_length_scale,false);
        }
        else
        {
            Ksat = fieldgen->interpolate("Ksat", actual_depth_outer);
            alpha = fieldgen->interpolate("alpha", actual_depth_outer);
            n = fieldgen->interpolate("n", actual_depth_outer);
            theta_s = fieldgen->interpolate("theta_s", actual_depth_outer);
            theta_r = fieldgen->interpolate("theta_r", actual_depth_outer);
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

            // NEW: Ksat scale factor
            B.SetVal("K_sat_scale_factor", KsatScale);

            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);

            // Geometry
            double actX = (i+0.5)*dr + mp.rw_uw;
            double actY = -actual_depth;

            B.SetVal("x",-(i*dr+mp.rw_uw)*2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",actX);
            B.SetVal("act_Y",actY);

            // INITIAL THETA
            double theta_init = theta_from_ert(actX, actual_depth, mp.initial_theta);
            B.SetVal("theta",theta_init);

            B.SetVal("L",mp.L);

            system->AddBlock(B,false);
        }

        // Soil Blocks directly under well
        std::cout<<"Soil Blocks directly under well"<<std::endl;

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

            // NEW: Ksat scale factor
            B.SetVal("K_sat_scale_factor", KsatScale);

            B.SetVal("area",area);
            B.SetVal("_width",dr*500);
            B.SetVal("_height",dr*500);
            B.SetVal("bottom_elevation",-(j+1)*dz-mp.DepthofWell_t);
            B.SetVal("depth",dz);

            // Geometry
            double actX = 0.0;
            double actY = -actual_depth;

            B.SetVal("x",-mp.rw_uw*1000+2000);
            B.SetVal("y",37000+(j*dz)*2000);
            B.SetVal("act_X",actX);
            B.SetVal("act_Y",actY);

            // INITIAL THETA
            double theta_init = theta_from_ert(actX, actual_depth, mp.initial_theta);
            B.SetVal("theta",theta_init);

            B.SetVal("L",mp.L);

            system->AddBlock(B,false);
        }
    }

    // --- remainder unchanged (links, wells, rain, solver setup) ---

    // Horizontal links for soils of gravel part
    std::cout<<"Horizontal links for soils of gravel part"<<std::endl;
    for (int i=0; i<mp.nr_g; i++)
        for (int j=0; j<mp.nz_g; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");

            L.SetName(("HL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ") - Soil-g (" + QString::number(i+2) + "$" + QString::number(j)+ ")").toStdString());
            L.SetType("soil_to_soil_H_link");

            L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_g));
            L.SetVal("length",dr);

            system->AddLink(L,
                ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(),
                ("Soil-g (" + QString::number(i+2) + "$" + QString::number(j) + ")").toStdString(),
                false);
        }

    // Vertical links for soils of gravel part
    std::cout<<"Vertical links for soils of gravel part"<<std::endl;
    for (int i=0; i<mp.nr_g; i++)
        for (int j=0; j<mp.nz_g-1; j++)
        {
            Link L;
            L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

            L.SetName(("VL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ") - Soil-g (" + QString::number(i+1) + "$" + QString::number(j+1)+ ")").toStdString());
            L.SetType("soil_to_soil_link");

            system->AddLink(L,
                ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(),
                ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j+1) + ")").toStdString(),
                false);
        }

    // Horizontal links for soils under well
    std::cout<<"Horizontal links for soils under well"<<std::endl;
    for (int i=0; i<mp.nr_uw; i++)
        for (int j = 0; j < nz_uw_n; j++)
            if (j*dz<mp.DepthtoGroundWater)
            {
                Link L;
                L.SetQuantities(system->GetMetaModel(), "soil_to_soil_H_link");

                L.SetName(("HL-Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ") - Soil-uw (" + QString::number(i+1) + "$" + QString::number(j)+ ")").toStdString());
                L.SetType("soil_to_soil_H_link");

                L.SetVal("area",2*pi*((i+0.5)*dr+mp.rw_uw));
                L.SetVal("length",dr);

                system->AddLink(L,
                    ("Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ")").toStdString(),
                    ("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j) + ")").toStdString(),
                    false);
            }

    // Vertical links for soils under well (first one)
    std::cout<<"Vertical links for soils under well (first one)"<<std::endl;
    for (int i=0; i<mp.nr_uw; i++)
    {
        Link L;
        L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

        int j_g=mp.nz_g-1;
        int j_uw=0;
        L.SetName(("VL-Soil-g (" + QString::number(i+1) + "$" + QString::number(j_g) + ") - Soil-uw (" + QString::number(i+1) + "$" + QString::number(j_uw)+ ")").toStdString());
        L.SetType("soil_to_soil_link");

        system->AddLink(L,
            ("Soil-g (" + QString::number(i+1) + "$" + QString::number(j_g) + ")").toStdString(),
            ("Soil-uw (" + QString::number(i+1) + "$" + QString::number(j_uw) + ")").toStdString(),
            false);
    }

    // Vertical links for soils under well
    std::cout<<"Vertical links for soils under well"<<std::endl;
    for (int i=0; i<mp.nr_uw+1; i++)
        for (int j = 0; j < nz_uw_n; j++)
            if (j*dz<mp.DepthtoGroundWater)
            {
                Link L;
                L.SetQuantities(system->GetMetaModel(), "soil_to_soil_link");

                L.SetName(("VL-Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ") - Soil-uw (" + QString::number(i) + "$" + QString::number(j+1)+ ")").toStdString());
                L.SetType("soil_to_soil_link");

                system->AddLink(L,
                    ("Soil-uw (" + QString::number(i) + "$" + QString::number(j) + ")").toStdString(),
                    ("Soil-uw (" + QString::number(i) + "$" + QString::number(j+1) + ")").toStdString(),
                    false);
            }

    // Well_c
    std::cout<<"Well_c"<<std::endl;
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
    if (raincfg.rain_data ==5)
        well_c.SetProperty("inflow", rain_file);
    system->AddBlock(well_c,false);

    // Well_g
    std::cout<<"Well_g"<<std::endl;
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

    // Well to well overflow
    std::cout<<"Well to well overflow"<<std::endl;
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

    // Junction
    std::cout<<"Junction"<<std::endl;
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

    // Well to junction
    std::cout<<"Well to junction"<<std::endl;
    Link well_to_junction;
    well_to_junction.SetQuantities(system->GetMetaModel(), "darcy_connector");
    well_to_junction.SetName("Well_to_junction");
    well_to_junction.SetType("darcy_connector");
    system->AddLink(well_to_junction,"Well_c","Junction_elastic",false);

    // Junction to well
    std::cout<<"Junction to well"<<std::endl;
    Link junction_to_well;
    junction_to_well.SetQuantities(system->GetMetaModel(), "darcy_connector");
    junction_to_well.SetName("Junction_to_well");
    junction_to_well.SetType("darcy_connector");
    system->AddLink(junction_to_well,"Junction_elastic","Well_g",false);

    // Well to gravel links
    std::cout<<"Horizontal links for well to gravel part soils"<<std::endl;

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

            system->AddLink(L, "Well_g",
                ("Soil-g (" + QString::number(i_g) + "$" + QString::number(j) + ")").toStdString(),
                false);
        }

    // Vertical links for well to under well soils
    std::cout<<"Vertical links for well to under well soils"<<std::endl;
    Link L1;
    L1.SetQuantities(system->GetMetaModel(), "Well2soil vertical link");

    int i_uw=0;
    int j_uw=0;
    L1.SetName(("VL_Well_g - Soil-uw (" + QString::number(i_uw) + "$" + QString::number(j_uw)+ ")").toStdString());
    L1.SetType("Well2soil vertical link");

    system->AddLink(L1, "Well_g",
        ("Soil-uw (" + QString::number(i_uw) + "$" + QString::number(j_uw) + ")").toStdString(),
        false);

    // Graoundwater fixed head
    std::cout<<"Groundwater"<<std::endl;
    Block gw;
    gw.SetQuantities(system->GetMetaModel(), "fixed_head");
    gw.SetName("Ground Water");
    gw.SetType("fixed_head");
    gw.SetVal("_height",500);
    gw.SetVal("_width",mp.RadiousOfInfluence*1000);
    gw.SetVal("head",-mp.DepthtoGroundWater);
    gw.SetVal("Storage",100000);
    gw.SetVal("x",-mp.nr_uw*1000);
    gw.SetVal("y",37000+(mp.nz_c+mp.nz_g+nz_uw_n)*2000);
    system->AddBlock(gw,false);

    // Groundwater links
    std::cout<<"Soil to Groundwater"<<std::endl;
    for (int i=0; i<mp.nr_uw+1; i++)
    {
        Link L2;
        L2.SetQuantities(system->GetMetaModel(), "soil_to_fixedhead_link");
        L2.SetName(("Soil to Groundwater (" + QString::number(i) + ")").toStdString());
        L2.SetType("soil_to_fixedhead_link");

        system->AddLink(L2,
            ("Soil-uw (" + QString::number(i) + "$" + QString::number(nz_uw_n-1) + ")").toStdString(),
            "Ground Water",
            false);
    }

    if (raincfg.rain_data != 5)
    {
        // Rain
        std::cout<<"Rain"<<std::endl;
        Source rain;
        rain.SetQuantities(system->GetMetaModel(), "Precipitation");
        rain.SetType("Precipitation");
        rain.SetName("Rain");
        rain.SetVal("_height",3000);
        rain.SetVal("_width",3000);
        rain.SetVal("x",-5000);
        rain.SetVal("y",-500);

        rain.SetProperty("timeseries", rain_file);
        system->AddSource(rain, false);

        // Catchment
        std::cout<<"Catchment"<<std::endl;
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
        std::cout<<"Catchment to well link"<<std::endl;

        // Catchment to well
        std::cout<<"Catchment to well"<<std::endl;
        Link L3;
        L3.SetQuantities(system->GetMetaModel(), "Surface water to well");
        L3.SetName("Catchment to well");
        L3.SetType("Surface water to well");
        L3.SetVal("ManningCoeff",mp.ManningCoeff_cmw);
        L3.SetVal("length",mp.length_cmw);
        system->AddLink(L3, "Catchment", "Well_c", false);
    }

    // Tracer (mean age)
    std::cout<<"Tracer"<<std::endl;
    Constituent mean_age_tracer;
    mean_age_tracer.SetQuantities(system->GetMetaModel(), "Constituent");
    mean_age_tracer.SetName("meanagetracer");
    mean_age_tracer.SetType("Constituent");
    system->AddConstituent(mean_age_tracer,false);

    // Reaction (aging)
    std::cout<<"Reaction"<<std::endl;
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

    system->SetSettingsParameter("maximum_time_allowed",simcfg.maximum_time_allowed);
    system->SetSettingsParameter("maximum_number_of_matrix_inverstions",simcfg.maximum_matrix_inversions);

    system->SetSystemSettings();

    std::cout<<"Populate functions"<<std::endl;
    system->PopulateOperatorsFunctions();
    std::cout<<"Variable parents"<<std::endl;
    system->SetVariableParents();
    return true;
}
