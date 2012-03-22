#include "sphere_surface.h"
#include "io_stl.h"
#include <cmath>

using namespace sloc;

static double A = 0.525731112119133606;
static double B = 0.850650808352039932;

static unsigned int octahedron_indices[8][3] = {
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 4},
    {0, 4, 1},
    {5, 2, 1},
    {5, 3, 2},
    {5, 4, 3},
    {5, 1, 4}
};

static double octahedron_vertices[6][3] = {
    { 0,  0, -1},
    { 1,  0,  0},
    { 0, -1,  0},
    {-1,  0,  0},
    { 0,  1,  0},
    { 0,  0,  1}
};

static unsigned int icosahedron_indices[20][3] = {
    {0, 4, 1},
    {0, 9, 4},
    {9, 5, 4},
    {4, 5, 8},
    {4, 8, 1},
    {8, 10, 1},
    {8, 3, 10},
    {5, 3, 8},
    {5, 2, 3},
    {2, 7, 3},
    {7, 10, 3},
    {7, 6, 10},
    {7, 11, 6},
    {11, 0, 6},
    {0, 1, 6},
    {6, 1, 10},
    {9, 0, 11},
    {9, 11, 2},
    {9, 2, 5},
    {7, 2, 11}
};

static double icosahedron_vertices[12][3] = {
    { A,  0, -B},
    {-A,  0, -B},
    { A,  0,  B},
    {-A,  0,  B},
    { 0, -B, -A},
    { 0, -B,  A},
    { 0,  B, -A},
    { 0,  B,  A},
    {-B, -A,  0},
    { B, -A,  0},
    {-B,  A,  0},
    { B,  A,  0}
};

inline static void normalize_vertex(double *a)
{
    double d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    a[0] /= d;
    a[1] /= d;
    a[2] /= d;
}


// ----------------------------------------------------------------------------

SphereSurface::SphereSurface() : radius(1.0) {}

SphereSurface::SphereSurface(double radius) : radius(radius) {}

SphereSurface::~SphereSurface() {}

void SphereSurface::set_radius(double radius)
{
    this->radius = radius;
}

void SphereSurface::set_refinement(unsigned int levels)
{
    set_refinement(levels, ICOSAHEDRON_BASE);
}

void SphereSurface::set_refinement(unsigned int levels, SphereSurface::polyhedron_base_t base)
{
    switch (base)
    {
        case OCTAHEDRON_BASE:
            make_sphere_with_octahedron_base(levels);
            break;
        case ICOSAHEDRON_BASE:
            make_sphere_with_icosahedron_base(levels);
            break;
    }
}

void SphereSurface::write_stl(const char *filename)
{
    unsigned int e;
    unsigned int nel = _facets.size() / 3;

    STL_File stl;
    stl.set_facets(nel);
    for (e = 0; e < nel; e++)
    {
        const unsigned int a = _facets[3*e+0];
        const unsigned int b = _facets[3*e+1];
        const unsigned int c = _facets[3*e+2];
        stl.set_facet_vertex(e, 0, _vertices[3*a+0], _vertices[3*a+1], _vertices[3*a+2]);
        stl.set_facet_vertex(e, 1, _vertices[3*b+0], _vertices[3*b+1], _vertices[3*b+2]);
        stl.set_facet_vertex(e, 2, _vertices[3*c+0], _vertices[3*c+1], _vertices[3*c+2]);
    }
    stl.write(filename);
}

// ----------------------------------------------------------------------------

void SphereSurface::make_sphere_with_octahedron_base(unsigned int levels)
{
    for (unsigned int i = 0; i < 8; i++)
    {
        subdivide_triangle(
            octahedron_vertices[octahedron_indices[i][0]],
            octahedron_vertices[octahedron_indices[i][1]],
            octahedron_vertices[octahedron_indices[i][2]],
            levels);
    }
}

void SphereSurface::make_sphere_with_icosahedron_base(unsigned int levels)
{
    for (unsigned int i = 0; i < 20; i++)
    {
        subdivide_triangle(
            icosahedron_vertices[icosahedron_indices[i][0]],
            icosahedron_vertices[icosahedron_indices[i][1]],
            icosahedron_vertices[icosahedron_indices[i][2]],
            levels);
    }
}

void SphereSurface::subdivide_triangle(double *a, double *b, double *c, unsigned int level)
{
    unsigned int i;
    double ab[3], ac[3], bc[3];

    if (level > 0)
    {
        //
        //          c
        //
        //      ac      bc
        //
        //  a       ab        b
        //

        for (i = 0; i < 3; i++)
        {
            ab[i] = (a[i] + b[i])/2;
            ac[i] = (a[i] + c[i])/2;
            bc[i] = (b[i] + c[i])/2;
        }

        normalize_vertex(ab);
        normalize_vertex(ac);
        normalize_vertex(bc);

        subdivide_triangle( a, ab, ac, level-1);
        subdivide_triangle( b, bc, ab, level-1);
        subdivide_triangle( c, ac, bc, level-1);
        subdivide_triangle(ab, bc, ac, level-1);
    }
    else
    {
        // add one more facet
        unsigned int n = _vertices.size() / 3;
        _facets.push_back(n);
        _facets.push_back(n+1);
        _facets.push_back(n+2);

        // add vertex n+1 with coords r*a and normal a
        for (i = 0; i < 3; i++) _normals.push_back(a[i]);
        for (i = 0; i < 3; i++) _vertices.push_back(radius * a[i]);

        // add vertex n+2 with coords r*b and normal b
        for (i = 0; i < 3; i++) _normals.push_back(b[i]);
        for (i = 0; i < 3; i++) _vertices.push_back(radius * b[i]);

        // add vertex n+3 with coords r*c and normal c
        for (i = 0; i < 3; i++) _normals.push_back(c[i]);
        for (i = 0; i < 3; i++) _vertices.push_back(radius * c[i]);
    }
}

