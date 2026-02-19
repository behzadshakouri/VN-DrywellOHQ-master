#ifndef MODELCREATOR_H
#define MODELCREATOR_H

#include <gsl/gsl_rng.h>

class System;
class FieldGenerator;

enum class _realization_mode { deterministic, stochastic };

// =============================================================
//   Simulation Control Structs (NEW)
// =============================================================

struct RainConfig
{
    int rain_data = 4;      // Default Pacoima; main() will set
};

struct SimulationConfig
{
    double start_time = 0.0;
    double end_time   = 1.0;

    double Base_start = 0.0;
    double Base_end   = 1.0;

    double maximum_time_allowed      = 10 * 86400;     // default 10 days
    double maximum_matrix_inversions = 10 * 200000;    // default 2,000,000
};

// =============================================================
//   Model Parameters (unchanged except removing limits)
// =============================================================

struct model_parameters
{
    double ft = 0.3048;
    double in = 0.0254;

    double DepthtoGroundWater = 142 * ft;
    double DepthofWell_c      = 16 * ft;
    double DepthofWell_g      = 24 * ft;
    double DepthofWell_t      = DepthofWell_c + DepthofWell_g;

    double RadiousOfInfluence = 20;

    int nr_c  = 4;
    int nr_g  = 12;
    int nr_uw = 12;

    int nz_c   = 5;
    int nz_g   = 15;
    int nz_uw  = 30;
    int nz_uw_n = 12; // For Synthetic rain (ERT_3: 13-12=1; ERT_5: 24-12=12)

    // Soil parameters
    double theta_r   = 0.049;
    double theta_sat = 0.39;
    double alpha     = 3.47536;
    double n         = 1.74582;
    double K_sat     = 1.05196;
    double K_o       = 0.24322;
    double L         = -0.5;

    double rw_c  = 4 * ft;
    double rw_c_t = 6 * ft;
    double rw_g  = 4 * ft;
    double rw_uw = 4 * ft;

    double initial_theta = 0.2;
    double porosity_c    = 1.0;
    double porosity_g    = 0.5;

    // Pipes
    double ManningCoeff_of   = 0.01;
    double diameter_of       = 8 * in;
    double length_of         = 10;
    double start_elevation_of = -6 * ft;
    double end_elevation_of   = -28 * ft;

    // Catchment
    double ManningCoeff_cm  = 0.03;
    double Slope_cm         = 0.01;
    double area_cm          = 2700;
    double Width_cm         = 15;
    double ManningCoeff_cmw = 0.05;
    double length_cmw       = 50;

    // Junction
    double elevation_j = -1 * DepthofWell_c;

    // Water depths
    double depth_w_c = 0.0;
    double depth_w_g = 0.01;

    double correlation_length_scale = 1.0;
};

// =============================================================
//   Initial-theta (ERT) mode switches (NEW)
// =============================================================
//
// These are intentionally small/simple ints so you can toggle
// them from main without pulling extra headers.
//
// ERTInitialThetaMode:
//   0 -> ERT-based through r (IDW in radius), depth-varying profile
//   1 -> ERT-based average through r (one depth-varying profile)
//
// ERTInitialThetaWhich:
//   0 -> use BOTH ERT-3 and ERT-5
//   3 -> use ONLY ERT-3
//   5 -> use ONLY ERT-5
//

// =============================================================
//   ModelCreator Class
// =============================================================

class ModelCreator
{
public:
    model_parameters ModelParameters() { return modelparameters; }

    ModelCreator();

    bool Create(model_parameters mp,
                System* system,
                FieldGenerator* fieldgen,
                const RainConfig& raincfg,
                const SimulationConfig& simcfg);

    // ------------------------------------------------------------
    // Optional initial-theta assignment from ERT tidy CSVs
    // (ERT-3_tidy.csv, ERT-5_tidy.csv). Uses the FIRST time slice
    // in each tidy file and maps theta [%] -> theta [0..1].
    // ------------------------------------------------------------
    bool UseERTInitialTheta = false;

    // NEW: Strategy switches (read by modelcreator.cpp)
    int ERTInitialThetaMode  = 0; // 0=IDW through r, 1=average through r
    int ERTInitialThetaWhich = 0; // 0=both, 3=ERT-3 only, 5=ERT-5 only

    _realization_mode Mode = _realization_mode::stochastic;

private:
    const double pi = 3.141521;
    model_parameters modelparameters;
};

#endif // MODELCREATOR_H
