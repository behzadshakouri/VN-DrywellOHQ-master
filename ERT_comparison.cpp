#include "ERT_comparison.h"

// std
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <cctype>

// Need full definition of model_parameters fields used below
#include "modelcreator.h"

#include "System.h"
#include "resultgrid.h"

// ------------------------------------------------------------
// Shared helpers (now shared with main.cpp via header prototypes)
// ------------------------------------------------------------
bool file_exists_cpp(const std::string& p)
{
    std::ifstream f(p);
    return (bool)f;
}

std::string lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

// ------------------------------------------------------------
// CLI parsing for init-theta
// ------------------------------------------------------------

// NEW: moved-from-main helper (declaration in ERT_comparison.h)
bool cli_has_init_theta(int argc, char** argv)
{
    // user explicitly provided init-theta?
    //   --init-theta <val>
    //   --init-theta=<val>
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i] ? argv[i] : "";
        if (a == "--init-theta") return true;
        if (a.rfind("--init-theta=", 0) == 0) return true;
    }
    return false;
}

InitThetaMode parse_init_theta_mode(int argc, char** argv)
{
    // Default should mean "use mp.initial_theta" unless user overrides via CLI.
    InitThetaMode mode = InitThetaMode::Default;

    auto normalize = [](std::string v)
    {
        v = lower_copy(v);
        // allow a couple separators
        std::replace(v.begin(), v.end(), '-', '_');
        return v;
    };

    auto apply_token = [&](const std::string& vraw)
    {
        const std::string v = normalize(vraw);

        if (v == "default" || v == "off" || v == "none" || v == "mp")
            mode = InitThetaMode::Default;
        else if (v == "ert3" || v == "ert3_only" || v == "3")
            mode = InitThetaMode::ERT3_Only;
        else if (v == "ert5" || v == "ert5_only" || v == "5")
            mode = InitThetaMode::ERT5_Only;
        else if (v == "idw" || v == "ert" || v == "r" || v == "ert_idw_r" || v == "idw_r")
            mode = InitThetaMode::ERT_IDW_R;
        else if (v == "avg" || v == "mean" || v == "ravg" || v == "ert_r_avg" || v == "r_avg")
            mode = InitThetaMode::ERT_R_Avg;
        // else: ignore unknown tokens (keep prior mode)
    };

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i] ? std::string(argv[i]) : std::string();

        const std::string keyEq = "--init-theta=";
        if (a.rfind(keyEq, 0) == 0)
        {
            apply_token(a.substr(keyEq.size()));
        }
        else if (a == "--init-theta" && i + 1 < argc)
        {
            apply_token(argv[++i]);
        }
        // Optional convenience flags (do not require a value)
        else if (a == "--ert3")     mode = InitThetaMode::ERT3_Only;
        else if (a == "--ert5")     mode = InitThetaMode::ERT5_Only;
        else if (a == "--ert-idw")  mode = InitThetaMode::ERT_IDW_R;
        else if (a == "--ert-ravg") mode = InitThetaMode::ERT_R_Avg;
    }

    return mode;
}

const char* init_theta_mode_name(InitThetaMode m)
{
    switch (m)
    {
        case InitThetaMode::Default:   return "Default";
        case InitThetaMode::ERT3_Only: return "ERT3_Only";
        case InitThetaMode::ERT5_Only: return "ERT5_Only";
        case InitThetaMode::ERT_IDW_R: return "ERT_IDW_R";
        case InitThetaMode::ERT_R_Avg: return "ERT_R_Avg";
        default: return "Unknown";
    }
}

// -----------------------------------------------------------------------------
// Legacy radius->index mapping (for legacy pipeline only)
// IMPORTANT: Soil-uw indices in your model are:
//   i=0      -> center column (act_X = 0)
//   i=1..nr  -> annular rings outside rw_uw
// -----------------------------------------------------------------------------
int mapRadiusToSoilUwIndex(double r_m, const model_parameters& mp)
{
    if (mp.nr_uw <= 0) return 0;
    if (!std::isfinite(r_m)) return 0;
    if (r_m <= 0.5 * mp.rw_uw) return 0;           // near center -> i=0
    if (r_m <= mp.rw_uw) return 1;                 // inside well radius -> closest annulus

    const double dr = (mp.RadiousOfInfluence - mp.rw_uw) / double(mp.nr_uw);
    if (dr <= 0) return 1;

    // annulus i=1 has "band" [rw_uw, rw_uw+dr), etc.
    int i = 1 + (int)std::floor((r_m - mp.rw_uw) / dr);
    i = std::clamp(i, 1, (int)mp.nr_uw);
    return i;
}

// -----------------------------------------------------------------------------
// NEW: use actual act_X rings from the model (Soil-uw / Soil-g)
// -----------------------------------------------------------------------------
struct RingInfo {
    int i;       // radial index in Soil-XX (i$j)
    double r;    // act_X radius (m)
};

static std::vector<RingInfo> get_rings_from_system(System* system,
                                                  const std::string& prefix,
                                                  int max_i)
{
    std::vector<RingInfo> rings;
    if (!system) return rings;

    // Use j=0 as representative (act_X should be constant per ring)
    for (int i = 0; i <= max_i; ++i)
    {
        std::string bname = prefix + " (" + std::to_string(i) + "$0)";
        auto* b = system->block(bname);
        if (!b) continue;

        double r = b->GetVal("act_X");
        if (!std::isfinite(r)) continue;

        rings.push_back({ i, r });
    }

    std::sort(rings.begin(), rings.end(),
              [](const RingInfo& a, const RingInfo& b){ return a.r < b.r; });

    return rings;
}

static int pick_nearest_ring_index(double r_m, const std::vector<RingInfo>& rings)
{
    if (rings.empty()) return 0;

    int best_i = rings.front().i;
    double best_d = std::abs(r_m - rings.front().r);

    for (const auto& rg : rings)
    {
        double d = std::abs(r_m - rg.r);
        if (d < best_d) { best_d = d; best_i = rg.i; }
    }
    return best_i;
}

static bool pick_bracketing_rings(double r_m, const std::vector<RingInfo>& rings,
                                 int& i0, int& i1, double& w)
{
    if (rings.empty()) {
        i0 = i1 = 0; w = 0.0; return false;
    }
    if (rings.size() == 1) {
        i0 = i1 = rings[0].i; w = 0.0; return true;
    }

    if (r_m <= rings.front().r) { i0 = i1 = rings.front().i; w = 0.0; return true; }
    if (r_m >= rings.back().r)  { i0 = i1 = rings.back().i;  w = 0.0; return true; }

    auto it = std::lower_bound(
        rings.begin(), rings.end(), r_m,
        [](const RingInfo& a, double x){ return a.r < x; }
    );

    const auto& hi = *it;
    const auto& lo = *(it - 1);

    i0 = lo.i;
    i1 = hi.i;

    double denom = (hi.r - lo.r);
    if (std::abs(denom) < 1e-30) { w = 0.0; return true; }

    w = (r_m - lo.r) / denom;
    w = std::clamp(w, 0.0, 1.0);
    return true;
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

    while (std::getline(f, line))
    {
        if (line.empty()) continue;

        if (first)
        {
            first = false;
            bool has_alpha = std::any_of(line.begin(), line.end(), [](unsigned char c){
                return std::isalpha(c);
            });
            if (has_alpha) continue; // skip header line
        }

        std::stringstream ss(line);
        std::string a, b;
        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;

        try {
            double d  = std::stod(a) + depth_offset_m;
            double th = std::stod(b);
            depth_m.push_back(d);
            theta_obs.push_back(th);
        } catch (...) {
            // skip bad lines
        }
    }

    // sort by depth
    if (depth_m.size() >= 2)
    {
        std::vector<size_t> idx(depth_m.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

        std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j){
            return depth_m[i] < depth_m[j];
        });

        std::vector<double> d2, t2;
        d2.reserve(depth_m.size());
        t2.reserve(theta_obs.size());

        for (size_t k = 0; k < idx.size(); ++k)
        {
            d2.push_back(depth_m[idx[k]]);
            t2.push_back(theta_obs[idx[k]]);
        }
        depth_m.swap(d2);
        theta_obs.swap(t2);
    }

    return !depth_m.empty();
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
    if (it == x.end())   return y.back();

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

    if (!bh.obs_csv_path.empty())
    {
        if (file_exists_cpp(bh.obs_csv_path))
            have_obs = read_observed_profile_csv(bh.obs_csv_path, d_obs, th_obs, bh.obs_depth_offset_m);

        if (!have_obs && !opt.allow_missing_obs)
            throw std::runtime_error("Observed CSV missing/unreadable: " + bh.obs_csv_path);
    }

    // Clamp time index using any series that exists
    unsigned int jt = opt.time_index;
    if (uniformoutput.size() > 0)
    {
        for (unsigned int k = 0; k < uniformoutput.size(); ++k)
        {
            if (uniformoutput[k].size() > 0)
            {
                if (jt >= uniformoutput[k].size())
                    jt = (unsigned int)(uniformoutput[k].size() - 1);
                break;
            }
        }
    }

    // Attempt to get time value from the first matching theta series if possible
    double t = 0.0;
    {
        std::string probeName = "Soil-uw (" + std::to_string(i_uw) + "$0)_" + opt.model_theta_var;
        for (unsigned int k = 0; k < uniformoutput.size(); ++k)
        {
            if (uniformoutput.getSeriesName(k) == probeName && uniformoutput[k].size() > 0)
            {
                if (jt >= uniformoutput[k].size())
                    jt = (unsigned int)(uniformoutput[k].size() - 1);
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
        std::string seriesName =
            "Soil-uw (" + std::to_string(i_uw) + "$" + std::to_string(j) + ")_" + opt.model_theta_var;

        double th_m = std::numeric_limits<double>::quiet_NaN();
        for (unsigned int k = 0; k < uniformoutput.size(); ++k)
        {
            if (uniformoutput.getSeriesName(k) == seriesName)
            {
                if (uniformoutput[k].size() > 0)
                {
                    unsigned int jtc = jt;
                    if (jtc >= uniformoutput[k].size())
                        jtc = (unsigned int)(uniformoutput[k].size() - 1);

                    th_m = uniformoutput[k].getValue(jtc);
                    t    = uniformoutput[k].getTime(jtc);
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

    // Rings from System for BOTH zones
    auto rings_g  = get_rings_from_system(system, "Soil-g",  mp.nr_g);
    auto rings_uw = get_rings_from_system(system, "Soil-uw", mp.nr_uw);

    if (rings_g.empty() && rings_uw.empty())
        throw std::runtime_error("No Soil-g / Soil-uw rings found in System (act_X missing?)");

    // Switch here if you ever want nearest-only
    const bool radial_interpolate = true;

    int i0_g=0, i1_g=0; double w_g=0.0;
    int i0_u=0, i1_u=0; double w_u=0.0;

    if (!rings_g.empty()) {
        if (radial_interpolate) pick_bracketing_rings(bh.r_m, rings_g,  i0_g, i1_g, w_g);
        else { i0_g = i1_g = pick_nearest_ring_index(bh.r_m, rings_g); w_g = 0.0; }
    }
    if (!rings_uw.empty()) {
        if (radial_interpolate) pick_bracketing_rings(bh.r_m, rings_uw, i0_u, i1_u, w_u);
        else { i0_u = i1_u = pick_nearest_ring_index(bh.r_m, rings_uw); w_u = 0.0; }
    }

    // Location lists for ALL blocks along depth at the selected rings
    std::vector<std::string> locs_g0, locs_g1;
    std::vector<std::string> locs_u0, locs_u1;

    if (!rings_g.empty()) {
        locs_g0.reserve((size_t)mp.nz_g);
        locs_g1.reserve((size_t)mp.nz_g);
        for (int j = 0; j < mp.nz_g; ++j) {
            locs_g0.push_back("Soil-g (" + std::to_string(i0_g) + "$" + std::to_string(j) + ")");
            locs_g1.push_back("Soil-g (" + std::to_string(i1_g) + "$" + std::to_string(j) + ")");
        }
    }

    if (!rings_uw.empty()) {
        locs_u0.reserve((size_t)opt.nz_uw_n);
        locs_u1.reserve((size_t)opt.nz_uw_n);
        for (int j = 0; j < opt.nz_uw_n; ++j) {
            locs_u0.push_back("Soil-uw (" + std::to_string(i0_u) + "$" + std::to_string(j) + ")");
            locs_u1.push_back("Soil-uw (" + std::to_string(i1_u) + "$" + std::to_string(j) + ")");
        }
    }

    // Build ResultGrids once (fast)
    ResultGrid rg_g0, rg_g1, rg_u0, rg_u1;

    if (!locs_g0.empty()) {
        rg_g0 = ResultGrid(uniformoutput, locs_g0, opt.model_theta_var);
        if (rg_g0.size() == 0 || rg_g0[0].size() == 0)
            throw std::runtime_error("rg_g0 empty/time-empty for borehole " + bh.name);

        if (i1_g != i0_g) {
            rg_g1 = ResultGrid(uniformoutput, locs_g1, opt.model_theta_var);
            if (rg_g1.size() == 0 || rg_g1[0].size() == 0) {
                rg_g1 = rg_g0; i1_g = i0_g; w_g = 0.0;
            }
        } else {
            rg_g1 = rg_g0;
        }
    }

    if (!locs_u0.empty()) {
        rg_u0 = ResultGrid(uniformoutput, locs_u0, opt.model_theta_var);
        if (rg_u0.size() == 0 || rg_u0[0].size() == 0)
            throw std::runtime_error("rg_u0 empty/time-empty for borehole " + bh.name);

        if (i1_u != i0_u) {
            rg_u1 = ResultGrid(uniformoutput, locs_u1, opt.model_theta_var);
            if (rg_u1.size() == 0 || rg_u1[0].size() == 0) {
                rg_u1 = rg_u0; i1_u = i0_u; w_u = 0.0;
            }
        } else {
            rg_u1 = rg_u0;
        }
    }

    // Clamp time index using whichever grid exists
    unsigned int jt = opt.time_index;
    if (!locs_g0.empty()) {
        if (jt >= rg_g0[0].size()) jt = (unsigned int)(rg_g0[0].size() - 1);
    } else {
        if (jt >= rg_u0[0].size()) jt = (unsigned int)(rg_u0[0].size() - 1);
    }

    // Time stamp
    double t = (!locs_g0.empty()) ? rg_g0[0].getTime(jt) : rg_u0[0].getTime(jt);

    // Read observed profile if provided
    std::vector<double> d_obs, th_obs;
    bool have_obs = false;

    if (!bh.obs_csv_path.empty())
    {
        if (file_exists_cpp(bh.obs_csv_path))
            have_obs = read_observed_profile_csv(bh.obs_csv_path, d_obs, th_obs, bh.obs_depth_offset_m);

        if (!have_obs && !opt.allow_missing_obs)
            throw std::runtime_error("Observed CSV missing/unreadable: " + bh.obs_csv_path);
    }

    // Collect rows then sort by depth
    struct Row {
        double depth;
        double th_m;
        double th_o;
    };
    std::vector<Row> rows;
    rows.reserve((size_t)mp.nz_g + (size_t)opt.nz_uw_n);

    auto obs_at_depth = [&](double depth)->double {
        double th_o = std::numeric_limits<double>::quiet_NaN();
        if (!have_obs) return th_o;
        if (d_obs.size() >= 2) return linear_interp_or_nan(d_obs, th_obs, depth);
        if (d_obs.size() == 1 && std::abs(d_obs[0] - depth) < 1e-6) return th_obs[0];
        return th_o;
    };

    // Collect Soil-g rows
    if (!locs_g0.empty())
    {
        for (int j = 0; j < mp.nz_g; ++j)
        {
            const std::string& bname = locs_g0[(size_t)j];
            auto* b = system->block(bname);
            if (!b) continue;

            double actY = b->GetVal("act_Y");
            if (!std::isfinite(actY)) continue;

            double depth = -actY; // act_Y is negative depth in your model

            double th0 = rg_g0[(size_t)j].getValue(jt);
            double th1 = rg_g1[(size_t)j].getValue(jt);

            double th_m = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(th0) && std::isfinite(th1)) th_m = (1.0 - w_g) * th0 + w_g * th1;
            else if (std::isfinite(th0)) th_m = th0;
            else if (std::isfinite(th1)) th_m = th1;

            rows.push_back({depth, th_m, obs_at_depth(depth)});
        }
    }

    // Collect Soil-uw rows
    if (!locs_u0.empty())
    {
        for (int j = 0; j < opt.nz_uw_n; ++j)
        {
            const std::string& bname = locs_u0[(size_t)j];
            auto* b = system->block(bname);
            if (!b) continue;

            double actY = b->GetVal("act_Y");
            if (!std::isfinite(actY)) continue;

            double depth = -actY;

            double th0 = rg_u0[(size_t)j].getValue(jt);
            double th1 = rg_u1[(size_t)j].getValue(jt);

            double th_m = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(th0) && std::isfinite(th1)) th_m = (1.0 - w_u) * th0 + w_u * th1;
            else if (std::isfinite(th0)) th_m = th0;
            else if (std::isfinite(th1)) th_m = th1;

            rows.push_back({depth, th_m, obs_at_depth(depth)});
        }
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b){ return a.depth < b.depth; });

    std::string outPath = opt.out_dir;
    if (!outPath.empty() && outPath.back() != '/' && outPath.back() != '\\') outPath += "/";
    outPath += opt.file_prefix + bh.name + opt.file_suffix;

    std::ofstream f(outPath);
    if (!f) throw std::runtime_error("Cannot open output: " + outPath);

    f << "depth_m,theta_model,theta_obs,time\n";
    f << std::setprecision(15);

    for (const auto& r : rows)
    {
        f << r.depth << "," << r.th_m << ",";
        if (std::isfinite(r.th_o)) f << r.th_o;
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

// ============================================================================
// moved-from-main: ERT snapshot driver + helpers
// ============================================================================

std::vector<double> default_ert_measurement_times()
{
    return {
        45763.588194444,
        45763.682638889,
        45763.770833333,
        45764.673611111,
        45764.729166667
    };
}

std::vector<BoreholeSpec> default_ert_boreholes(const std::string& workingFolder)
{
    std::string wf = workingFolder;
    if (!wf.empty() && wf.back() != '/' && wf.back() != '\\') wf += "/";

    return {
        {"ERT-3",  6.7979, wf + "obs/ERT-3_obs.csv", 0.0},
        {"ERT-5",  4.3974, wf + "obs/ERT-5_obs.csv", 0.0},
    };
}

std::string ert_time_token(double t_excel_days)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << t_excel_days;
    std::string s = oss.str();
    std::replace(s.begin(), s.end(), '.', 'p');
    return s;
}

unsigned int nearest_time_index_in_uniform_set(
    const TimeSeriesSet<double>& uniformoutput,
    double target_t_excel_days)
{
    if (uniformoutput.size() == 0) return 0;

    const auto& ts0 = uniformoutput[0];
    unsigned int n = (unsigned int)ts0.size();
    if (n == 0) return 0;

    unsigned int best_i = 0;
    double best_dt = std::numeric_limits<double>::infinity();

    for (unsigned int i = 0; i < n; ++i)
    {
        double t  = ts0.getTime(i);
        double dt = std::abs(t - target_t_excel_days);
        if (dt < best_dt) { best_dt = dt; best_i = i; }
    }
    return best_i;
}

void export_ert_snapshots_at_measurement_times(
    const TimeSeriesSet<double>& uniformoutput_ERT,
    const model_parameters& mp,
    System* system,
    const RainConfig& raincfg,
    bool use_resultgrid_ert,
    bool allow_missing_obs)
{
    if (!system) throw std::runtime_error("System is null");

    const std::string wf = system->GetWorkingFolder();
    std::vector<BoreholeSpec> boreholes = default_ert_boreholes(wf);
    std::vector<double> ert_times = default_ert_measurement_times();

    // same nz selection you used in main
    int nz_uw_n = (raincfg.rain_data == 5) ? mp.nz_uw_n : mp.nz_uw;

    BoreholeExportOptions opt;
    opt.nz_uw_n           = nz_uw_n;
    opt.out_dir           = wf;
    opt.model_theta_var   = "theta";
    opt.allow_missing_obs = allow_missing_obs;
    opt.use_resultgrid    = use_resultgrid_ert;

    for (double t_ert : ert_times)
    {
        unsigned int idx = nearest_time_index_in_uniform_set(uniformoutput_ERT, t_ert);
        std::string tag  = "t" + ert_time_token(t_ert);

        opt.file_prefix = "ERTsnap-" + tag + "-";
        opt.file_suffix = ".csv";
        opt.time_index  = idx;

        std::cout << "\n=== ERT snapshot export ===\n";
        std::cout << "Target ERT time = " << std::setprecision(12) << t_ert << "\n";
        std::cout << "Nearest model index (ERT grid) = " << idx << "\n";
        std::cout << "Output prefix = " << opt.file_prefix << "\n";

        export_borehole_theta_comparison_csvs(
            uniformoutput_ERT, mp, system,
            boreholes, opt
        );
    }
}
