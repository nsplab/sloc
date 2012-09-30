#ifndef SLOC_SPATIAL_INDEX_H
#define SLOC_SPATIAL_INDEX_H

#include <vector>
//#include <set>

#include <boost/shared_ptr.hpp>
#include <octree/point3d.h>
#include <octree/octree.h>
#include "bbox3.h"

namespace sloc
{

    template <int N = 1024>
    class SpatialIndex
    {
    public:

        typedef nomis80::Point3D<double> point_t;
        typedef nomis80::Point3D<int> voxel_t;
        typedef std::vector<int> bucket_t;
        typedef boost::shared_ptr<bucket_t> bucket_ptr;

    public:

        SpatialIndex(sloc::bbox3& bbox);
        ~SpatialIndex();

        bool find_voxel(const point_t& pt, voxel_t& voxel);
        void add_points(const std::vector<point_t>& xs);
        int add_point(point_t p);
        int next_index() const;
        int find_index(const point_t& p);

    public:

        sloc::bbox3& bbox;
        std::vector<point_t> points;
        nomis80::Octree<bucket_ptr> octree;

    };


} // namespace sloc

#endif
