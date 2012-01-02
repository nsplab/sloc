#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <vector>

namespace sloc
{
// ----------------------------------------------------------------------------

class PointCloud
{
public:
    PointCloud(double tol);
    ~PointCloud();
    void clear();

    int n_points() const;
    void get_point(int n, double *point);
    void add(double x, double y, double z, long *id);

public:
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
