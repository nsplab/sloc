#include "sloc/io_getfem.h"
#include "sloc/io_stl.h"
#include <fstream>

using namespace std;

void sloc::read_stl(getfem::mesh& mesh, const char *filename)
{
    //
    // Given the data in the stl file, we need to create
    // a mesh object with all duplicate points removed.
    //
    // Luckily, getfem::mesh removes the duplicate points, so
    // we just load the points into that structure.
    //
    const bool verbose = false;

    if (verbose)
        cout << "sloc::stl_read() filename=" << filename << endl;

    STL_File stl;
    stl.read(filename);

    int ncells = stl.n_facets();
    int npts = ncells * 3;

    if (verbose)
    {
        cout << "  stl facets = " << ncells << endl;
        cout << "  stl points = " << npts << endl;
    }

    // loop indices
    int e, n;

    // loop over STL facets
    for (e = 0; e < ncells; e++)
    {
        vector<bgeot::size_type> ind(3);
        ind[0] = mesh.add_point(bgeot::base_node(stl.va[3*e+0], stl.va[3*e+1], stl.va[3*e+2]));
        ind[1] = mesh.add_point(bgeot::base_node(stl.vb[3*e+0], stl.vb[3*e+1], stl.vb[3*e+2]));
        ind[2] = mesh.add_point(bgeot::base_node(stl.vc[3*e+0], stl.vc[3*e+1], stl.vc[3*e+2]));
        mesh.add_convex(bgeot::simplex_geotrans(2,1), ind.begin());
    }

}

void sloc::read_points(std::vector<bgeot::base_node>& points, const char *filename)
{
    std::ifstream is;
    is.open(filename);

    unsigned int n = 0;
    is >> n;

    for (unsigned int i = 0; i < n; ++i)
    {
        double x, y, z;
        is >> x >> y >> z;
        points.push_back(bgeot::base_node(x,y,z));
    }

    is.close();
}


/*
void sloc::write_points(const std::vector<bgeot::base_node>& points, const char *filename)
{
    std::ofstream os;
    os.open(filename);

    const unsigned int n = points.size();
    os << n << endl;
    for (unsigned int i = 0; i < n; ++i)
    {
        const bgeot::base_node& pt = points[i];
        os << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    }
}
*/

