#ifndef RESULTGRID_H
#define RESULTGRID_H

#include "TimeSeries.h"
#include "TimeSeriesSet.h"

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
    ResultGrid(const TimeSeriesSet<double> &cts, const vector<string> &components, const string &quantity);

    TimeSeries<double> Sum();
    TimeSeries<double> SumIntegrate();
    vector<point> Positions;
    void WriteToVTP(const std::string &quanname, const std::string &filename, int i, const double &scale=1) const;
    void WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale=1) const;

};

#endif // RESULTGRID_H
