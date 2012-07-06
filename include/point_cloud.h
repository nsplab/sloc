#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <vector>
#include <octree/octree.h>
#include <boost/smart_ptr.hpp>

namespace sloc {

class PointCloud
{
public:

    PointCloud();
    ~PointCloud();

    void clear();

    inline int n_points() const { return npts; }

    void set_tolerance(double tol);

    void create_index(int size);

    void add(double x, double y, double z, long *id);

    void get_point(int n, double *point);

    bool naive_search(double x, double y, double z, long *id);
    bool search(double x, double y, double z, long *id);

public:

    double epsilon;

    int npts;
    std::vector<double*> points;

    typedef std::vector<int> bucket_t;
    typedef boost::shared_ptr<bucket_t> shared_bucket_t;
    typedef nomis80::Octree<shared_bucket_t> octree_idx;
    boost::shared_ptr<octree_idx> idx;
};

} // namespace sloc
#endif
