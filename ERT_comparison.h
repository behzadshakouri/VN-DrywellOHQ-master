#pragma once

#include <string>
#include <vector>

// Forward declarations (avoid heavy includes in header)
class System;
template <typename T> class TimeSeriesSet;

struct model_parameters;
struct RainConfig;

// ============================================================
// Small shared helpers (used by main.cpp too)
// ============================================================
bool file_exists_cpp(const std::string& p);
std::string lower_copy(std::string s);

// ============================================================
// CLI: Initial-theta mode (used by main.cpp to configure ModelCreator)
// ============================================================
enum class InitThetaMode
{
    Default = 0,
    ERT3_Only,
    ERT5_Only,
    ERT_IDW_R,
    ERT_R_Avg
};

InitThetaMode parse_init_theta_mode(int argc, char** argv);
const char* init_theta_mode_name(InitThetaMode m);

// ============================================================
// Data structures for ERT/Borehole exports
// ============================================================

struct BoreholeSpec
{
    std::string name;                 // e.g., "ERT-3"
    double r_m = 0.0;                 // radial distance from well centerline (m)
    std::string obs_csv_path;         // optional CSV path for observations
    double obs_depth_offset_m = 0.0;  // optional offset added to observed depths
};

struct BoreholeExportOptions
{
    int nz_uw_n = 0;
    unsigned int time_index = 0;
    std::string out_dir;
    std::string file_prefix = "BH_";
    std::string file_suffix = "_theta_compare.csv";
    std::string model_theta_var = "theta";
    bool allow_missing_obs = true;
    bool use_resultgrid = true;
};

// ============================================================
// Core export
// ============================================================

void export_borehole_theta_comparison_csvs(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const std::vector<BoreholeSpec>& boreholes,
    const BoreholeExportOptions& opt);

// ============================================================
// Higher-level "ERT snapshots @ measurement times" driver
// ============================================================

std::vector<double> default_ert_measurement_times();
std::vector<BoreholeSpec> default_ert_boreholes(const std::string& workingFolder);

std::string ert_time_token(double t_excel_days);

unsigned int nearest_time_index_in_uniform_set(
    const TimeSeriesSet<double>& uniformoutput,
    double target_t_excel_days);

void export_ert_snapshots_at_measurement_times(
    const TimeSeriesSet<double>& uniformoutput_ERT,
    const model_parameters& mp,
    System* system,
    const RainConfig& raincfg,
    bool use_resultgrid_ert = true,
    bool allow_missing_obs  = true);

// ============================================================
// Utilities
// ============================================================

int mapRadiusToSoilUwIndex(double r_m, const model_parameters& mp);

bool read_observed_profile_csv(
    const std::string& path,
    std::vector<double>& depth_m,
    std::vector<double>& theta_obs,
    double depth_offset_m = 0.0);

double linear_interp_or_nan(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double xq);
