#ifndef RESULTGRID_H
#define RESULTGRID_H

#include "TimeSeries.h"
#include "TimeSeriesSet.h"

// ----------------------------------------------------------------------------
// ResultGrid
//
// Lightweight spatial wrapper over TimeSeriesSet<double>.
//
// Concept:
//   - Each TimeSeries<double> represents a variable at one spatial location.
//   - ResultGrid adds spatial awareness by maintaining a parallel
//     Positions[] vector aligned with the internal TimeSeriesSet indices.
//
// Important:
//   Positions[i] corresponds to the spatial location of
//   (*this)[i]  (i.e., the i-th TimeSeries).
//
// Supports optional VTK export when compiled with use_VTK.
// ----------------------------------------------------------------------------

// Project VTK include (your build system likely redirects this)
#include <vtk.h>

#include <string>
#include <vector>
#include <map>

// ----------------------------------------------------------------------------
// Simple 2D spatial point
// ----------------------------------------------------------------------------
struct point
{
    double x = 0.0;   // X coordinate
    double y = 0.0;   // Y coordinate
};

class System;

// ----------------------------------------------------------------------------
// ResultGrid
// ----------------------------------------------------------------------------
// Inherits from TimeSeriesSet<double> and augments it with spatial metadata.
// Used for:
//   - Storing distributed model outputs
//   - Snapshot extraction at specific times
//   - CSV export (ERT / Borehole comparisons)
//   - Optional 3D VTK visualization
// ----------------------------------------------------------------------------
class ResultGrid : public TimeSeriesSet<double>
{
public:

    // ------------------------------------------------------------------------
    // Constructors / Rule of Three
    // ------------------------------------------------------------------------

    ResultGrid();                              // Default constructor
    ResultGrid(const ResultGrid&);              // Copy constructor
    virtual ~ResultGrid();                      // Destructor
    ResultGrid& operator=(const ResultGrid&);   // Copy assignment

    // Construct from an existing TimeSeriesSet and attach quantity name + system
    ResultGrid(const TimeSeriesSet<double> &cts,
               const std::string &quantity,
               System *system);

    // Construct empty grid with quantity label
    ResultGrid(const std::string &quantity,
               System *system);

    // Construct with spatial filtering (x-range)
    ResultGrid(const TimeSeriesSet<double> &cts,
               const string &quantity,
               const double &x_min,
               const double &x_max,
               System *system);

    // Construct multi-component quantity (e.g., vector fields)
    ResultGrid(const TimeSeriesSet<double> &cts,
               const std::vector<std::string> &components,
               const std::string &quantity);

    // ------------------------------------------------------------------------
    // Aggregate utilities
    // ------------------------------------------------------------------------

    // Sum across all spatial series (pointwise in time)
    // Returns one aggregated TimeSeries
    TimeSeries<double> Sum();

    // Time-integrated sum across all spatial series
    // Typically used for cumulative mass/energy quantities
    TimeSeries<double> SumIntegrate();

    // ------------------------------------------------------------------------
    // Spatial metadata
    // ------------------------------------------------------------------------

    // Spatial coordinates aligned with TimeSeries indices.
    // Must always remain synchronized with inherited container.
    std::vector<point> Positions;

    // ------------------------------------------------------------------------
    // Snapshots & CSV export (ERT / Borehole helpers)
    // ------------------------------------------------------------------------

    // Snapshot at a single time t (interpolated).
    //
    // Returns:
    //   A new ResultGrid where each spatial cell contains
    //   a single-point TimeSeries evaluated at time t.
    ResultGrid Snapshot(const double &t) const;

    // Snapshot at multiple times (interpolated).
    //
    // Returns:
    //   A new ResultGrid where each series is sampled at the
    //   provided time vector.
    ResultGrid Snapshot(const std::vector<double> &t) const;

    // Write snapshot values to CSV using a given time index.
    //
    // CSV columns:
    //   name, x, y, <value_col>, time
    //
    // time_index:
    //   Index into the internal time vector of each series.
    //
    // value_col:
    //   Column name for exported quantity (default = "value").
    bool WriteSnapshotCSV(const std::string& filename,
                          unsigned int time_index,
                          const std::string& value_col = "value") const;

#ifdef use_VTK

    // ------------------------------------------------------------------------
    // VTK EXPORT FUNCTIONS
    // ------------------------------------------------------------------------

    // Write a single snapshot to VTP file.
    //
    // i      : time index
    // scale  : scaling factor applied to geometry (visual scaling)
    void WriteToVTPSnapShot(const std::string &quanname,
                            const std::string &filename,
                            int i,
                            const double &scale) const;

    // Write full time sequence to VTP files (animation-ready)
    //
    // starting_counter:
    //   Offset for file numbering
    void WriteToVTP(const std::string &quanname,
                    const std::string &filename,
                    const double &scale,
                    int starting_counter) const;

    // ------------------------------------------------------------------------
    // Geometry Builders (Static Helpers)
    // ------------------------------------------------------------------------

    // Create a solid cylinder
    static vtkSmartPointer<vtkPolyData>
    MakeCylinder(double radius,
                 double height,
                 double centerZ);

    // Create hollow cylinder (outer/inner radii)
    // values: map of scalar attributes attached to geometry
    static vtkSmartPointer<vtkPolyData>
    MakeHollowCylinder(double zCenter,
                       double height,
                       double outerRadius,
                       double innerRadius,
                       std::map<std::string, float> values);

    // Build full 3D VTK dataset from multiple quantities
    static void Make3DVTK(std::vector<std::string> quantity,
                          double dr,
                          System *system,
                          std::string filename);

    // Solid tube
    static vtkSmartPointer<vtkPolyData>
    MakeTube(double radius,
             double height,
             double zCenter,
             std::map<std::string, float> values);

    // Hollow tube (outer & inner radii)
    static vtkSmartPointer<vtkPolyData>
    MakeTube(double outerRadius,
             double innerRadius,
             double height,
             double zCenter,
             std::map<std::string, float> values);

    // Manual hollow tube builder (custom implementation)
    static vtkSmartPointer<vtkPolyData>
    MakeHollowTubeManual(double outerRadius,
                         double innerRadius,
                         double height,
                         double zCenter,
                         std::map<std::string, float> values);

    // Alternative hollow tube builder
    static vtkSmartPointer<vtkPolyData>
    MakeHollowTube(double outerRadius,
                   double innerRadius,
                   double height,
                   double zCenter,
                   std::map<std::string, float> values);

#endif // use_VTK
};

#endif // RESULTGRID_H
