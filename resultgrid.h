#ifndef RESULTGRID_H
#define RESULTGRID_H

#include "TimeSeries.h"
#include "TimeSeriesSet.h"
#include <vtk.h>

struct point
{
    double x;
    double y;
};

class System;

class ResultGrid: public TimeSeriesSet<double>
{
public:
    ResultGrid();
    ResultGrid(const ResultGrid&);
    virtual ~ResultGrid();
    ResultGrid& operator=(const ResultGrid&);
    ResultGrid(const TimeSeriesSet<double> &cts, const string &quantity, System *system);
    ResultGrid(const string &quantity, System *system);
    ResultGrid(const TimeSeriesSet<double> &cts, const vector<string> &components, const string &quantity);

    TimeSeries<double> Sum();
    TimeSeries<double> SumIntegrate();
    vector<point> Positions;
    void WriteToVTP(const std::string &quanname, const std::string &filename, int i, const double &scale=1) const;
    void WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale=1) const;

    static vtkSmartPointer<vtkPolyData> MakeCylinder(double radius, double height, double centerZ);
    static vtkSmartPointer<vtkPolyData> MakeHollowCylinder(double zCenter, double height, double outerRadius, double innerRadius,map<string, float> values);
    static void Make3DVTK(vector<string> quantity, double dr,System *system, string filename);
    static vtkSmartPointer<vtkPolyData> MakeTube(double radius, double height, double zCenter, std::map<std::string, float> values);
    static vtkSmartPointer<vtkPolyData> MakeTube(double outerRadius, double innerRadius, double height, double zCenter, std::map<std::string, float> values);
    static vtkSmartPointer<vtkPolyData> MakeHollowTubeManual(double outerRadius, double innerRadius, double height, double zCenter, std::map<std::string, float> values);
    static vtkSmartPointer<vtkPolyData> MakeHollowTube(double outerRadius, double innerRadius, double height, double zCenter, std::map<std::string, float> values);

};

#endif // RESULTGRID_H
