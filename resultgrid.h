#ifndef RESULTGRID_H
#define RESULTGRID_H

#include "TimeSeries.h"
#include "TimeSeriesSet.h"

// ----------------------------------------------------------------------------
// ResultGrid
// - Lightweight spatial wrapper over TimeSeriesSet<double>
// - Keeps a parallel Positions[] vector aligned with series indices
// - Supports VTK export when compiled with use_VTK
// ----------------------------------------------------------------------------

// Keep this include as you had it (your project likely routes VTK includes here)
#include <vtk.h>

#include <string>
#include <vector>
#include <map>

struct point
{
    double x = 0.0;
    double y = 0.0;
};

class System;

class ResultGrid : public TimeSeriesSet<double>
{
public:
    ResultGrid();
    ResultGrid(const ResultGrid&);
    virtual ~ResultGrid();
    ResultGrid& operator=(const ResultGrid&);

    ResultGrid(const TimeSeriesSet<double> &cts, const std::string &quantity, System *system);
    ResultGrid(const std::string &quantity, System *system);
    ResultGrid(const TimeSeriesSet<double> &cts, const std::vector<std::string> &components, const std::string &quantity);

    // Aggregate utilities
    TimeSeries<double> Sum();
    TimeSeries<double> SumIntegrate();

    // Spatial positions (must stay aligned with series indices)
    std::vector<point> Positions;

    // ------------------------------------------------------------------------
    // NEW: Snapshots and simple CSV export (ERT helpers)
    // ------------------------------------------------------------------------

    // Returns a ResultGrid with 1-point TimeSeries per cell at time t (interpolated).
    ResultGrid Snapshot(const double &t) const;

    // Returns a ResultGrid where each series is sampled at given times (interpolated).
    ResultGrid Snapshot(const std::vector<double> &t) const;

    // Writes one CSV of a snapshot at a specific index (per-series value at that time index).
    // CSV columns: name,x,y,<value_col>,time
    bool WriteSnapshotCSV(const std::string& filename,
                          unsigned int time_index,
                          const std::string& value_col = "value") const;

#ifdef use_VTK
    // VTK exports
    void WriteToVTPSnapShot(const std::string &quanname, const std::string &filename, int i, const double &scale) const;
    void WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale, int starting_counter) const;

    // Geometry builders
    static vtkSmartPointer<vtkPolyData> MakeCylinder(double radius, double height, double centerZ);

    // NOTE: This is declared in your old header. Keep it if used elsewhere.
    // If not implemented in .cpp, either implement or remove (but you said no deletions).
    static vtkSmartPointer<vtkPolyData> MakeHollowCylinder(double zCenter, double height,
                                                           double outerRadius, double innerRadius,
                                                           std::map<std::string, float> values);

    static void Make3DVTK(std::vector<std::string> quantity, double dr, System *system, std::string filename);

    static vtkSmartPointer<vtkPolyData> MakeTube(double radius, double height, double zCenter,
                                                 std::map<std::string, float> values);

    static vtkSmartPointer<vtkPolyData> MakeTube(double outerRadius, double innerRadius, double height, double zCenter,
                                                 std::map<std::string, float> values);

    static vtkSmartPointer<vtkPolyData> MakeHollowTubeManual(double outerRadius, double innerRadius, double height, double zCenter,
                                                             std::map<std::string, float> values);

    static vtkSmartPointer<vtkPolyData> MakeHollowTube(double outerRadius, double innerRadius, double height, double zCenter,
                                                       std::map<std::string, float> values);
#endif
};

#endif // RESULTGRID_H
