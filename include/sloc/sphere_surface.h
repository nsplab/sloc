#ifndef SLOC_SPHERE_SURFACE_H
#define SLOC_SPHERE_SURFACE_H

#include <vector>

namespace sloc {
// ----------------------------------------------------------------------------

class SphereSurface
{
public:

    typedef enum { ICOSAHEDRON_BASE, OCTAHEDRON_BASE } polyhedron_base_t;

    SphereSurface();
    SphereSurface(double r);
    ~SphereSurface();

    void set_radius(double radius);
    void set_refinement(unsigned int levels);
    void set_refinement(unsigned int levels, polyhedron_base_t base);
    void write_stl(const char *filename);

protected:

    void make_sphere_with_octahedron_base(unsigned int levels);
    void make_sphere_with_icosahedron_base(unsigned int levels);
    void subdivide_triangle(double *a, double *b, double *c, unsigned int level);

private:

    double radius;
    std::vector<double> _normals;
    std::vector<double> _vertices;
    std::vector<unsigned int> _facets;
};


// ----------------------------------------------------------------------------
}
#endif
