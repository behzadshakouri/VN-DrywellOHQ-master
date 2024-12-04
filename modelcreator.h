#ifndef MODELCREATOR_H
#define MODELCREATOR_H

#include <gsl/gsl_rng.h>

class System;

struct model_parameters
{
    double DepthtoGroundWater = 20; // 40'
    double DepthofWell = 12.1920; // 16'
    double DepthofWell_c = 4.8768; // 16'
    double DepthofWell_g1 = 0.6096; // 2'
    double DepthofWell_g2 = 3.6576; // 12'
    double RadiousOfInfluence = 1.2192; //4'
    int nr = 4; // Radius discretization number
    int nz = 10; // Depth discretization number
    double K_sat =1;
    double alpha = 20;
    double n = 1.8;
    double rw = 0.1;
    double rw_c = 1.2192; //4'
    double rw_g1 = 1.8288; //6'
    double rw_g2 = 1.2192; //4'
    double theta_sat = 0.4;
    double theta_r = 0.05;
    double initial_theta = 0.2;
    double porosity_c = 1;
    double porosity_g1 = 0.5;
    double porosity_g2 = 0.5;


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
