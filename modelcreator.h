#ifndef MODELCREATOR_H
#define MODELCREATOR_H

#include <gsl/gsl_rng.h>

class System;

struct model_parameters
{
    double ft = 0.3048; //ft to m
    double in = 0.0254; //in to m
    double DepthtoGroundWater = 85; // Should be estimated
    double DepthofWell_c = 16*ft; // 16'
    double DepthofWell_g = 24*ft; // 24'
    double DepthofWell_t = DepthofWell_c + DepthofWell_g; // 40'
    double RadiousOfInfluence = 20; //Should be estimated

    int nr_c = 4; // Radius discretization number (around concrete part)
    int nr_g = 12; // Radius discretization number (around gravel part)
    int nr_uw = 12; // Radius discretization number (under well part)

    int nz_c = 5; // Depth discretization number (around concrete part)
    int nz_g = 15; // Depth discretization number (around gravel part)
    int nz_uw = 30; // Depth discretization number (under well part)

    double K_sat = 0.075; //will be calculated
    double alpha = 5; //will be calculated
    double n = 1.8;
    double L = -0.874;
    double rw_c = 4*ft; //4'
    double rw_c_t = 6*ft; //6'
    double rw_g = 4*ft; //4'
    double rw_uw = 4*ft; //4'
    double theta_sat = 0.39;
    double theta_r = 0.049;
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

};

class ModelCreator
{
public:
    ModelCreator();
    bool Create(model_parameters mp, System *system);

private:
    const double pi = 3.141521;
};

#endif // MODELCREATOR_H
