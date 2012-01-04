#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <vector>

namespace sloc
{
// ----------------------------------------------------------------------------

class PointCloud
{
public:
    PointCloud();
    ~PointCloud();

    void clear();

    int n_points() const;

    void set_tolerance(double tol);

    void add(double x, double y, double z, long *id);

    void get_point(int n, double *point);

    bool naive_search(double x, double y, double z, long *id);

public:
    double epsilon;
    int npts;
    std::vector<double*> points;
};

inline int PointCloud::n_points() const { return npts; }

// ----------------------------------------------------------------------------
}
#endif
