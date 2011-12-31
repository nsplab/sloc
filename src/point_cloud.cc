#include "point_cloud.h"

using namespace std;

PointCloud::PointCloud(double tol)
{
    epsilon = tol;
    npts = 0;
}

PointCloud::~PointCloud()
{
    clear();
}

void PointCloud::clear()
{
    npts = 0;
    vector<double*>::iterator it;
    for (it = points.begin(); it != points.end(); ++it)
        delete *it;
}

void PointCloud::get_point(int n, double *point)
{
    // load n-th point
    double *p = points[n];

    // copy the coordinates
    point[0] = p[0];
    point[1] = p[1];
    point[2] = p[2];
}

void PointCloud::add(double x, double y, double z, long *id)
{
    bool found;

    // look for (x,y,z) in current point cloud
    found = naive_search(x, y, z, id);

    if (!found)
    {
        // new point!
        double *p = new double[3];
        p[0] = x; p[1] = y; p[2] = z;
        points.push_back(p);

        // use last index (zero-based) for the new point
        if (id != 0)
            *id = npts;

        // update count, since we've added a new point
        npts++;
    }
}

static double dist2(double x[3], double y[3])
{
    double d[3];
    d[0] = x[0] - y[0];
    d[1] = x[1] - y[1];
    d[2] = x[2] - y[2];
    return d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
}

bool PointCloud::naive_search(double x, double y, double z, long *id)
{
    double tol = epsilon * epsilon;

    double p[3];
    p[0] = x; p[1] = y; p[2] = z;

    for (unsigned long i = 0; i < points.size(); i++)
    {
        if (dist2(p, points[i]) < tol)
        {
            if (id != 0)
                *id = i;
            return true;
        }
    }

    return false;
}

// EOF
