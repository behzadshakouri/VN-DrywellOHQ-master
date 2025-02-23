#ifndef MODELCREATOR_H
#define MODELCREATOR_H

#include <gsl/gsl_rng.h>

class System;

struct model_parameters
{
    double ft=0.3048; //ft to m
    double DepthtoGroundWater = 20; // Should be estimated
    double DepthofWell_c = 16*ft; // 16'
    double DepthofWell_g = 24*ft; // 24'
    double DepthofWell_t = DepthofWell_c + DepthofWell_g; // 40'
    double RadiousOfInfluence = 20; //Should be estimated
    int nr_c = 4; // Radius discretization number (around concrete part)
    int nr_g = 8; // Radius discretization number (around gravel part)
    int nz_c = 5; // Depth discretization number (around concrete part)
    int nz_g = 10; // Depth discretization number (around gravel part)
    double K_sat =1;
    double alpha = 20;
    double n = 1.8;
    double rw_c_t = 6*ft; //6'
    double rw_c = 4*ft; //4'
    double rw_g = 4*ft; //4'
    double theta_sat = 0.4;
    double theta_r = 0.05;
    double initial_theta = 0.2;
    double porosity_c = 1;
    double porosity_g = 0.5;


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
