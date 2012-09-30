#include "spatial_index.h"

using namespace sloc;

template <int N>
SpatialIndex<N>::SpatialIndex(bbox3& bbox)
    : bbox(bbox)
    , octree(N, bucket_ptr())
{
    cout << "Initialized " << N << "x" << N << "x" << N << " octree with bounding box " << bbox << endl;
}

template <int N>
SpatialIndex<N>::~SpatialIndex()
{
}

template <int N>
bool SpatialIndex<N>::find_voxel(const point_t& pt, voxel_t& voxel)
{
    if (bbox.outside(pt(0), pt(1), pt(2)))
    {
        voxel(0) = -1;
        voxel(1) = -1;
        voxel(2) = -1;
        return false;
    }

    double x = pt(0);
    double y = pt(1);
    double z = pt(2);

    voxel(0) = static_cast<int>(floor(N * (x - bbox.min[0]) / (bbox.max[0] - bbox.min[0])));
    voxel(1) = static_cast<int>(floor(N * (y - bbox.min[1]) / (bbox.max[1] - bbox.min[1])));
    voxel(2) = static_cast<int>(floor(N * (z - bbox.min[2]) / (bbox.max[2] - bbox.min[2])));

    if (voxel(0) == N) { voxel(0)--; }
    if (voxel(1) == N) { voxel(1)--; }
    if (voxel(2) == N) { voxel(2)--; }

    return true;
}

template <int N>
void SpatialIndex<N>::add_points(const std::vector<point_t>& xs)
{
    vector<point_t>::const_iterator it;
    for (it = xs.begin(); it != xs.end(); ++it)
    {
        this->add_point(*it);
    }
}

template <int N>
void SpatialIndex<N>::add_point(point_t p)
{
    // which index?
    int idx = -1;

    voxel_t voxel;

    if (find_voxel(p, voxel))
    {
        // voxel components
        int i,j,k;
        i = voxel(0);
        j = voxel(1);
        k = voxel(2);

        // get bucket at voxel
        bucket_t *bucket = octree.at(i,j,k).get();

        // make one if it doesn't exist
        if (!bucket)
        {
            bucket = new bucket_t;
            octree(i,j,k) = bucket_ptr(bucket);
        }

        // if we just insert p into bucket we are allowing duplicates...
        // thus we have to search the current bucket and reuse that index
        for (bucket_t::iterator it = bucket->begin(); it != bucket->end(); ++it)
        {
            point_t& q = points[*it];
            if (close_enough(p,q))
            {
                idx = *it;
                break;
            }
        }

        if (idx == -1)
        {
            // we didn't find the point... add it to the bucket
            idx = next_index();
            bucket->push_back(idx);
            points.push_back(p);
        }
    
    }
    else
    {
        cerr << "error: point " << p << " lies outside bounding box " << bbox << endl;
    }

    return idx;
}

template <int N>
int SpatialIndex<N>::next_index() const
{
    // next available index is at the end of points list... so just return the size
    return points.size();
}

template <int N>
int SpatialIndex<N>::find_index(const point_t& p)
{
    // first, find the voxel
    voxel_t voxel;

    if (find_voxel(p, voxel))
    {
        // next, search the voxel's bucket
        bucket_t *bucket = octree.at(voxel(0), voxel(1), voxel(2)).get();

        if (bucket)
        {
            for (bucket_t::iterator it = bucket->begin(); it != bucket->end(); ++it)
            {
                if (close_enough(p, points[*it])) { return *it; }
            }
        }
    }
    return -1;
}

