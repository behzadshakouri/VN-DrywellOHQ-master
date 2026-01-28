#pragma once

#include <string>
#include <vector>
#include <stdexcept>

// Forward declarations (avoid heavy includes in header)
class System;
template <typename T> class TimeSeriesSet;
class ResultGrid;

// Your project's model_parameters type
// If it lives in another header, include that instead.
#include "modelcreator.h"   // for model_parameters in your codebase (adjust if different)

struct BoreholeSpec
{
    std::string name;            // e.g., "BH01"
    double r_m = 0.0;            // radial distance from well centerline (m)

    // Optional observed profile file (CSV columns: depth_m,theta_obs)
    // Leave empty if none.
    std::string obs_csv_path;

    // Optional: if your observed depths are relative to ground surface but model depths
    // are absolute, you can offset obs depths here (meters). Typically 0.
    double obs_depth_offset_m = 0.0;
};

struct BoreholeExportOptions
{
    // Under-well vertical layers count used during model creation:
    // nz_uw_n = (raincfg.rain_data == 5) ? mp.nz_uw_n : mp.nz_uw
    int nz_uw_n = 0;

    // Time sample index in uniformoutput (will be clamped to available range)
    unsigned int time_index = 0;

    // Output directory (ensure trailing slash if you want)
    std::string out_dir;

    // Output filename pattern:
    // out_dir + "BH_<name>_theta_compare.csv"
    std::string file_prefix = "BH_";
    std::string file_suffix = "_theta_compare.csv";

    // Model variable name for moisture
    std::string model_theta_var = "theta";

    // If true, missing obs file is not an error; theta_obs will be blank.
    bool allow_missing_obs = true;
};

// ------------------------------------------------------------
// Main API
// ------------------------------------------------------------

// Exports one CSV per borehole for theta(depth) comparison.
// CSV columns: depth_m, theta_model, theta_obs, time
void export_borehole_theta_comparison_csvs(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const std::vector<BoreholeSpec>& boreholes,
    const BoreholeExportOptions& opt);

// ------------------------------------------------------------
// Utilities (exposed in case you want to reuse)
// ------------------------------------------------------------

// Map borehole radius r (m) to Soil-uw radial index i (0..nr_uw)
int mapRadiusToSoilUwIndex(double r_m, const model_parameters& mp);

// Reads observed profile CSV depth_m,theta_obs.
// Returns true if successfully read; false if file missing/unreadable.
bool read_observed_profile_csv(
    const std::string& path,
    std::vector<double>& depth_m,
    std::vector<double>& theta_obs,
    double depth_offset_m = 0.0);

// Piecewise-linear interpolation (x must be sorted ascending).
// Returns NaN if xq outside range.
double linear_interp_or_nan(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double xq);
