/*
 * Create a sphere.surf.ucd file
 */

#include "sphere_surface.h"
#include "io_stl.h"
#include <iostream>

int main(void)
{
    using namespace std;
    using namespace sloc;

    const char *filename = "sphere.surf.stl";

    SphereSurface sphere_surf;
    sphere_surf.set_radius(1.0);
    sphere_surf.set_refinement(2);
    sphere_surf.write_stl(filename);

    STL_File stl;
    stl.read(filename);

    cout << "STL facets = " << stl.n_facets() << endl;
    if (stl.n_facets() < 100)
    {
        float x, y, z;
        int i,j;
        for (i = 0; i < stl.n_facets(); i++)
        {
            cout << i << " ";
            for (j = 0; j < 3; j++)
            {
                stl.get_facet_vertex(i, j, x, y, z);
                cout << "[" << x << ", " << y << ", " << z << "] ";
            }
            cout << endl;
        }
    }

    return 0;
}
