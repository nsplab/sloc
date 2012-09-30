#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <vector>
#include <octree/octree.h>
#include <octree/point3d.h>
#include <boost/smart_ptr.hpp>
#include "spatial_index.h"

namespace sloc
{

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
        bool octree_search(double x, double y, double z, long *id);

    public:

        double epsilon;

        int npts;
        std::vector<double*> points;

        typedef sloc::SpatialIndex<128> spatial_index_t;
        boost::shared_ptr<spatial_index_t> spatial_index;

    };

} // namespace sloc
#endif
