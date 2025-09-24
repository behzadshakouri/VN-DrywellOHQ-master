#ifndef MODELCREATOR_H
#define MODELCREATOR_H

#include <gsl/gsl_rng.h>

class System;
class FieldGenerator;

enum class _realization_mode {deterministic, stochastic};

struct model_parameters
{
    double ft = 0.3048; //ft to m
    double in = 0.0254; //in to m
    double DepthtoGroundWater = 85; // Should be estimated, 25 to 85 m
    double DepthofWell_c = 16*ft; // 16'
    double DepthofWell_g = 24*ft; // 24'
    double DepthofWell_t = DepthofWell_c + DepthofWell_g; // 40'
    double RadiousOfInfluence = 20; //Should be estimated

/*
    // optimized
    int nr_c = 4; // Radius discretization number (around concrete part), 4
    int nr_g = 12; // Radius discretization number (around gravel part), 12
    int nr_uw = 12; // Radius discretization number (under well part), 12

    int nz_c = 5; // Depth discretization number (around concrete part), 5
    int nz_g = 15; // Depth discretization number (around gravel part), 15
    int nz_uw = 30; // Depth discretization number (under well part), 30
*/


    // smaller for testing
    int nr_c = 4; // Radius discretization number (around concrete part), 4
    int nr_g = 6; // Radius discretization number (around gravel part), 12
    int nr_uw = 6; // Radius discretization number (under well part), 12

    int nz_c = 5; // Depth discretization number (around concrete part), 5
    int nz_g = 5; // Depth discretization number (around gravel part), 15
    int nz_uw = 8; // Depth discretization number (under well part), 30


    // Soil properties
    double theta_r = 0.049; // will be calculated, 0.049, Rosetta Sandy Loam Non-log 0.049
    double theta_sat = 0.39;// will be calculated, 0.39, Rosetta Sandy Loam Non-log 0.39
    double alpha = 3.47536; // will be calculated, 5, Rosetta Sandy Loam Non-log 3.47536
    double n = 1.74582; // will be calculated, 1.8, Rosetta Sandy Loam Non-log 1.74582
    double K_sat = 1.05196; // will be calculated, 0.3, Rosetta Sandy Loam Non-log 1.05196
    double K_o = 0.24322; // will be calculated, 0.24322, Rosetta Sandy Loam Non-log 0.24322
    double L = -0.5; // will be calculated, -0.874, Rosetta Sandy Loam Non-log -0.874, for new approach -0.5 (default)
    double rw_c = 4*ft; //4'
    double rw_c_t = 6*ft; //6'
    double rw_g = 4*ft; //4'
    double rw_uw = 4*ft; //4'
    double initial_theta = 0.2;
    double porosity_c = 1;
    double porosity_g = 0.5;
    double ManningCoeff_of = 0.01; //overflow
    double diameter_of = 8*in;
    double length_of = 10;
    double start_elevation_of = -2;
    double end_elevation_of = -10;
    double ManningCoeff_cm = 0.03; //catchment
    double Slope_cm = 0.01;
    double area_cm = 2700;
    double Width_cm = 15;
    double ManningCoeff_cmw = 0.05; //catchment to well link
    double length_cmw = 50;
    double elevation_j = -1*DepthofWell_c;
    double depth_w_c=0; //depth of water in well
    double depth_w_g=0.01; //depth of water in well
    double Maximum_time_allowed=10*86400; // 10 days
    double Maximum_number_of_matrix_inverstions=10*200000; // 10x

    double correlation_length_scale = 1.0;
};

class ModelCreator
{
public:
    model_parameters ModelParameters() {return modelparameters;}
    ModelCreator();
    bool Create(model_parameters mp, System *system, FieldGenerator *propgen = nullptr);
    _realization_mode Mode = _realization_mode::stochastic;
private:
    const double pi = 3.141521;
    model_parameters modelparameters;
};

#endif // MODELCREATOR_H
