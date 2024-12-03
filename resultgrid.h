#ifndef RESULTGRID_H
#define RESULTGRID_H

#include <BTCSet.h>

struct point
{
    double x;
    double y;
};

class System;

class ResultGrid: public CTimeSeriesSet<double>
{
public:
    ResultGrid();
    ResultGrid(const ResultGrid&);
    virtual ~ResultGrid();
    ResultGrid& operator=(const ResultGrid&);
    ResultGrid(const CTimeSeriesSet<double> &cts, const string &quantity, System *system);
    ResultGrid(const CTimeSeriesSet<double> &cts, const vector<string> &components, const string &quantity);
    CTimeSeries<double> Sum();
    CTimeSeries<double> SumIntegrate();
    vector<point> Positions;
    void WriteToVTP(const std::string &quanname, const std::string &filename, int i, const double &scale=1) const;
    void WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale=1) const;

};

#endif // RESULTGRID_H
