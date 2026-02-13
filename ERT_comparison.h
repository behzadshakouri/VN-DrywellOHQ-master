#pragma once

#include <string>
#include <vector>
#include <stdexcept>

// Forward declarations (avoid heavy includes in header)
class System;
template <typename T> class TimeSeriesSet;

// model_parameters lives in modelcreator.h in your codebase
#include "modelcreator.h"

/*
 * BoreholeSpec
 * ------------
 * Defines a borehole/ERT radial location in the model and an optional observed profile CSV.
 */
struct BoreholeSpec
{
    std::string name;                 // e.g., "ERT-2"
    double r_m = 0.0;                 // radial distance from well centerline (m)
    std::string obs_csv_path;         // optional CSV: depth_m,theta_obs
    double obs_depth_offset_m = 0.0;  // optional offset added to observed depths
};

/*
 * BoreholeExportOptions
 * ---------------------
 * Controls how profile exports are produced.
 */
struct BoreholeExportOptions
{
    int nz_uw_n = 0;                      // number of under-well vertical layers to export
    unsigned int time_index = 0;          // time index in uniformoutput (clamped)
    std::string out_dir;                  // output directory
    std::string file_prefix = "BH_";      // filename prefix
    std::string file_suffix = "_theta_compare.csv"; // filename suffix
    std::string model_theta_var = "theta";// variable name (quantity)
    bool allow_missing_obs = true;        // if obs file missing -> no error, theta_obs blank

    // NEW: switch between two internal pipelines (keeps older one intact).
    // true  -> uses ResultGrid-based extraction (recommended)
    // false -> uses legacy scanning from TimeSeriesSet (still works; kept as-is)
    bool use_resultgrid = true;
};

// Exports one CSV per borehole.
// CSV columns: depth_m,theta_model,theta_obs,time
void export_borehole_theta_comparison_csvs(
    const TimeSeriesSet<double>& uniformoutput,
    const model_parameters& mp,
    System* system,
    const std::vector<BoreholeSpec>& boreholes,
    const BoreholeExportOptions& opt);

// Utilities
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
