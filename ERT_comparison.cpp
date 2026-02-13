#include "ERT_comparison.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

#include "System.h"
#include "resultgrid.h"

static inline bool file_exists(const std::string& p)
{
    std::ifstream f(p);
    return (bool)f;
}

int mapRadiusToSoilUwIndex(double r_m, const model_parameters& mp)
{
    if (r_m <= 0.0) return 0;
    if (mp.nr_uw <= 0) return 0;

    const double dr = (mp.RadiousOfInfluence - mp.rw_uw) / double(mp.nr_uw);

    if (r_m <= mp.rw_uw) return 1;

    int i = 1 + (int)std::floor((r_m - mp.rw_uw) / dr);
    i = std::clamp(i, 1, (int)mp.nr_uw);
    return i;
}

bool read_observed_profile_csv(
    const std::string& path,
    std::vector<double>& depth_m,
    std::vector<double>& theta_obs,
    double depth_offset_m)
{
    depth_m.clear();
    theta_obs.clear();

    std::ifstream f(path);
    if (!f) return false;

    std::string line;
    bool first = true;

    while (std::getline(f, line)) {
        if (line.empty()) continue;

        if (first) {
            first = false;
            bool has_alpha = std::any_of(line.begin(), line.end(), [](unsigned char c){
                return std::isalpha(c);
            });
            if (has_alpha) continue;
        }

        std::stringstream ss(line);
        std::string a, b;
        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;

        try {
            double d = std::stod(a) + depth_offset_m;
            double th = std::stod(b);
            depth_m.push_back(d);
            theta_obs.push_back(th);
        } catch (...) {
            // skip bad lines
        }
    }

    if (depth_m.size() >= 2) {
        std::vector<size_t> idx(depth_m.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

        std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j){
            return depth_m[i] < depth_m[j];
        });

        std::vector<double> d2, t2;
        d2.reserve(depth_m.size());
        t2.reserve(theta_obs.size());

        for (size_t k = 0; k < idx.size(); ++k) {
            d2.push_back(depth_m[idx[k]]);
            t2.push_back(theta_obs[idx[k]]);
        }
        depth_m.swap(d2);
        theta_obs.swap(t2);
    }

    return true;
}

double linear_interp_or_nan(
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

// Legacy pipeline (kept inside this file; no change in outputs format)
static void export_one_borehole_legacy(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const BoreholeSpec& bh,
    const BoreholeExportOptions& opt)
{
    if (!system) throw std::runtime_error("System is null");

    const int i_uw = mapRadiusToSoilUwIndex(bh.r_m, mp);

    const double dz = (mp.DepthtoGroundWater - mp.DepthofWell_t) / double(mp.nz_uw);

    // Read observed profile if provided
    std::vector<double> d_obs, th_obs;
    bool have_obs = false;
    if (!bh.obs_csv_path.empty()) {
        have_obs = read_observed_profile_csv(bh.obs_csv_path, d_obs, th_obs, bh.obs_depth_offset_m);
        if (!have_obs && !opt.allow_missing_obs) {
            throw std::runtime_error("Observed CSV missing/unreadable: " + bh.obs_csv_path);
        }
    }

    // Clamp time index using any series that exists
    unsigned int jt = opt.time_index;
    if (uniformoutput.size() > 0) {
        // find first matching series for time support (best effort)
        for (unsigned int k = 0; k < uniformoutput.size(); ++k) {
            if (uniformoutput[k].size() > 0) {
                if (jt >= uniformoutput[k].size()) jt = (unsigned int)(uniformoutput[k].size() - 1);
                break;
            }
        }
    }

    // Attempt to get time value from the first matching theta series if possible
    double t = 0.0;
    {
        std::string probeName = "Soil-uw (" + std::to_string(i_uw) + "$0)_" + opt.model_theta_var;
        for (unsigned int k = 0; k < uniformoutput.size(); ++k) {
            if (uniformoutput.getSeriesName(k) == probeName && uniformoutput[k].size() > 0) {
                if (jt >= uniformoutput[k].size()) jt = (unsigned int)(uniformoutput[k].size() - 1);
                t = uniformoutput[k].getTime(jt);
                break;
            }
        }
    }

    std::string outPath = opt.out_dir;
    if (!outPath.empty() && outPath.back() != '/' && outPath.back() != '\\') outPath += "/";
    outPath += opt.file_prefix + bh.name + opt.file_suffix;

    std::ofstream f(outPath);
    if (!f) throw std::runtime_error("Cannot open output: " + outPath);

    f << "depth_m,theta_model,theta_obs,time\n";
    f << std::setprecision(15);

    for (int j = 0; j < opt.nz_uw_n; ++j)
    {
        const double depth = (j + 0.5) * dz + mp.DepthofWell_t;

        // Find the matching series by name
        std::string seriesName = "Soil-uw (" + std::to_string(i_uw) + "$" + std::to_string(j) + ")_" + opt.model_theta_var;

        double th_m = std::numeric_limits<double>::quiet_NaN();
        for (unsigned int k = 0; k < uniformoutput.size(); ++k) {
            if (uniformoutput.getSeriesName(k) == seriesName) {
                if (uniformoutput[k].size() > 0) {
                    unsigned int jtc = jt;
                    if (jtc >= uniformoutput[k].size()) jtc = (unsigned int)(uniformoutput[k].size() - 1);
                    th_m = uniformoutput[k].getValue(jtc);
                    t = uniformoutput[k].getTime(jtc);
                }
                break;
            }
        }

        double th_o = std::numeric_limits<double>::quiet_NaN();
        if (have_obs && d_obs.size() >= 2) th_o = linear_interp_or_nan(d_obs, th_obs, depth);
        else if (have_obs && d_obs.size() == 1 && std::abs(d_obs[0] - depth) < 1e-6) th_o = th_obs[0];

        f << depth << "," << th_m << ",";
        if (std::isfinite(th_o)) f << th_o;
        f << "," << t << "\n";
    }
}

// ResultGrid pipeline
static void export_one_borehole_resultgrid(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const BoreholeSpec& bh,
    const BoreholeExportOptions& opt)
{
    if (!system) throw std::runtime_error("System is null");
    if (opt.nz_uw_n <= 0) throw std::runtime_error("nz_uw_n <= 0");

    const double dz = (mp.DepthtoGroundWater - mp.DepthofWell_t) / double(mp.nz_uw);

    const int i_uw = mapRadiusToSoilUwIndex(bh.r_m, mp);

    std::vector<std::string> locs;
    locs.reserve((size_t)opt.nz_uw_n);
    for (int j = 0; j < opt.nz_uw_n; ++j) {
        locs.push_back("Soil-uw (" + std::to_string(i_uw) + "$" + std::to_string(j) + ")");
    }

    // Build a ResultGrid that contains theta time series for each depth location
    ResultGrid thetaRG(uniformoutput, locs, opt.model_theta_var);

    if (thetaRG.size() == 0) throw std::runtime_error("thetaRG empty for borehole " + bh.name);
    if (thetaRG[0].size() == 0) throw std::runtime_error("thetaRG time series empty for borehole " + bh.name);

    unsigned int jt = opt.time_index;
    if (jt >= thetaRG[0].size()) jt = (unsigned int)(thetaRG[0].size() - 1);

    const double t = thetaRG[0].getTime(jt);

    // Read observed profile if provided
    std::vector<double> d_obs, th_obs;
    bool have_obs = false;
    if (!bh.obs_csv_path.empty()) {
        have_obs = read_observed_profile_csv(bh.obs_csv_path, d_obs, th_obs, bh.obs_depth_offset_m);
        if (!have_obs && !opt.allow_missing_obs) {
            throw std::runtime_error("Observed CSV missing/unreadable: " + bh.obs_csv_path);
        }
    }

    std::string outPath = opt.out_dir;
    if (!outPath.empty() && outPath.back() != '/' && outPath.back() != '\\') outPath += "/";
    outPath += opt.file_prefix + bh.name + opt.file_suffix;

    std::ofstream f(outPath);
    if (!f) throw std::runtime_error("Cannot open output: " + outPath);

    f << "depth_m,theta_model,theta_obs,time\n";
    f << std::setprecision(15);

    for (int j = 0; j < opt.nz_uw_n; ++j)
    {
        const double depth = (j + 0.5) * dz + mp.DepthofWell_t;

        double th_m = thetaRG[j].getValue(jt);

        double th_o = std::numeric_limits<double>::quiet_NaN();
        if (have_obs && d_obs.size() >= 2) th_o = linear_interp_or_nan(d_obs, th_obs, depth);
        else if (have_obs && d_obs.size() == 1 && std::abs(d_obs[0] - depth) < 1e-6) th_o = th_obs[0];

        f << depth << "," << th_m << ",";
        if (std::isfinite(th_o)) f << th_o;
        f << "," << t << "\n";
    }
}

void export_borehole_theta_comparison_csvs(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const std::vector<BoreholeSpec>& boreholes,
    const BoreholeExportOptions& opt)
{
    if (!system) throw std::runtime_error("System is null");
    if (opt.nz_uw_n <= 0) throw std::runtime_error("nz_uw_n <= 0");
    if (opt.out_dir.empty()) throw std::runtime_error("out_dir empty");

    for (const auto& bh : boreholes)
    {
        if (opt.use_resultgrid) {
            export_one_borehole_resultgrid(uniformoutput, mp, system, bh, opt);
        } else {
            export_one_borehole_legacy(uniformoutput, mp, system, bh, opt);
        }
    }
}
