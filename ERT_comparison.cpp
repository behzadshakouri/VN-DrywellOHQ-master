#include "ERT_comparison.h""

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

// Your project includes (adjust paths if needed)
#include "System.h"
#include "resultgrid.h"

// -----------------------------
// Helpers
// -----------------------------

static inline bool file_exists(const std::string& p)
{
    std::ifstream f(p);
    return (bool)f;
}

int mapRadiusToSoilUwIndex(double r_m, const model_parameters& mp)
{
    // Valid indices in your model: i = 0..nr_uw (inclusive)
    // i=0 is the center column ("Soil-uw (0$j)")
    if (r_m <= 0.0) return 0;

    if (mp.nr_uw <= 0) return 0;
    const double dr = (mp.RadiousOfInfluence - mp.rw_uw) / double(mp.nr_uw);

    // For r <= rw_uw, closest ring is i=1 (first annulus outside the well)
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

        // skip header if present
        if (first) {
            first = false;
            // If header contains any alpha, treat as header line
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
            // ignore bad lines
        }
    }

    if (depth_m.size() < 2) return true; // still ok, but interp may be limited

    // sort by depth
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

// -----------------------------
// Core export
// -----------------------------

void export_borehole_theta_comparison_csvs(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const std::vector<BoreholeSpec>& boreholes,
    const BoreholeExportOptions& opt)
{
    if (!system) throw std::runtime_error("export_borehole_theta_comparison_csvs: system is null");
    if (opt.nz_uw_n <= 0) throw std::runtime_error("export_borehole_theta_comparison_csvs: opt.nz_uw_n <= 0");
    if (opt.out_dir.empty()) throw std::runtime_error("export_borehole_theta_comparison_csvs: opt.out_dir empty");

    // ModelCreator::Create under-well dz (same formula)
    const double dz = (mp.DepthtoGroundWater - mp.DepthofWell_t) / double(mp.nz_uw);

    for (const auto& bh : boreholes)
    {
        // Map radial distance to i index used in Soil-uw naming
        const int i_uw = mapRadiusToSoilUwIndex(bh.r_m, mp);

        // Build location list Soil-uw (i$j) for j=0..nz_uw_n-1
        std::vector<std::string> locs;
        locs.reserve((size_t)opt.nz_uw_n);
        for (int j = 0; j < opt.nz_uw_n; ++j) {
            locs.push_back("Soil-uw (" + std::to_string(i_uw) + "$" + std::to_string(j) + ")");
        }

        ResultGrid thetaRG(uniformoutput, locs, opt.model_theta_var);

        if (thetaRG.size() == 0)
            throw std::runtime_error("thetaRG empty for borehole " + bh.name +
                                     " (check naming or variable '" + opt.model_theta_var + "')");
        if (thetaRG[0].size() == 0)
            throw std::runtime_error("thetaRG time series empty for borehole " + bh.name);

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

        // Output CSV path
        std::string outPath = opt.out_dir;
        // Ensure separator
        if (!outPath.empty() && outPath.back() != '/' && outPath.back() != '\\')
            outPath += "/";

        outPath += opt.file_prefix + bh.name + opt.file_suffix;

        std::ofstream f(outPath);
        if (!f) throw std::runtime_error("Cannot open output: " + outPath);

        f << "depth_m,theta_model,theta_obs,time\n";
        f << std::setprecision(15);

        for (int j = 0; j < opt.nz_uw_n; ++j)
        {
            // depth used for those blocks:
            // actual_depth = (j+0.5)*dz + mp.DepthofWell_t
            const double depth = (j + 0.5) * dz + mp.DepthofWell_t;
            const double th_m = thetaRG[j].getValue(jt);

            // interpolate obs to model depths (common for comparison)
            double th_o = std::numeric_limits<double>::quiet_NaN();
            if (have_obs && d_obs.size() >= 2) {
                th_o = linear_interp_or_nan(d_obs, th_obs, depth);
            } else if (have_obs && d_obs.size() == 1) {
                // single-point obs: only match if very close
                if (std::abs(d_obs[0] - depth) < 1e-6) th_o = th_obs[0];
            }

            f << depth << "," << th_m << ",";
            if (std::isfinite(th_o)) f << th_o;
            f << "," << t << "\n";
        }
    }
}
