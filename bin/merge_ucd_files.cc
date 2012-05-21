/* merge_ucd_files.cc
 *
 * Use this program to merge Bone2_1.quad.ucd and Artery2_1.quad.ucd into head.ucd
 */

#include <iostream>
#include <string>
#include <cassert>
#include "mesh.h"

using namespace std;

// ----------------------------------------------------------------------------

void union_merge(sloc::Mesh& out, const sloc::Mesh& a, const sloc::Mesh& b)
{
    //
    // Assume nothing overlaps (no common nodes), and that our
    // two meshes have the same element type.
    //
    // We just copy the data over, and renumber the node ids.
    //

    cout << "Calling union_merge()" << endl;

    int npts = a.n_points() + b.n_points();
    int ncells = a.n_cells() + b.n_cells();

    int ncellnodes = a.n_cell_nodes();
    assert(a.n_cell_nodes() == b.n_cell_nodes());

    if (false)
    {
        cout << "  a.n_points()     = " << a.n_points() << endl;
        cout << "  a.n_cells()      = " << a.n_cells() << endl;
        cout << "  a.n_cell_nodes() = " << a.n_cell_nodes() << endl;

        cout << "  b.n_points()     = " << b.n_points() << endl;
        cout << "  b.n_cells()      = " << b.n_cells() << endl;
        cout << "  b.n_cell_nodes() = " << b.n_cell_nodes() << endl;

        cout << "  output npts      = " << npts << endl;
        cout << "  output ncells    = " << ncells << endl;
    }

    int e,n,i;
    int mat_id;
    double pt[3];
    long *cell = new long[a.n_cell_nodes()];

    out.init_points(npts, 3);
    out.init_cells(ncells, ncellnodes);

    //
    // copy a's data
    //

    for (n = 0; n < a.n_points(); n++)
    {
        a.get_point(n, pt);
        out.set_point(n, pt);
    }

    for (e = 0; e < a.n_cells(); e++)
    {
        a.get_cell(e, cell);
        out.set_cell(e, cell);
    }

    for (e = 0; e < a.n_cells(); e++)
    {
        a.get_mat(e, mat_id);
        out.set_mat(e, mat_id);
    }

    //
    // copy b's data, but renumber the node and cell ids
    //

    const int node_offset = a.n_points();
    const int cell_offset = a.n_cells();

    for (n = 0; n < b.n_points(); n++)
    {
        b.get_point(n, pt);
        out.set_point(n + node_offset, pt);
    }

    for (e = 0; e < b.n_cells(); e++)
    {
        b.get_cell(e, cell);

        for (i = 0; i < ncellnodes; i++)
            cell[i] += node_offset;

        out.set_cell(e + cell_offset, cell);
    }

    for (e = 0; e < b.n_cells(); e++)
    {
        b.get_mat(e, mat_id);
        out.set_mat(e + cell_offset, mat_id);
    }

    delete [] cell;

}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    sloc::Mesh merged_mesh, first_mesh, second_mesh;
    string merged_filename, first_filename, second_filename;

    if (argc == 1)
    {
        cout << "Usage: " << argv[0]
             << " <first-mesh.ucd> <second-mesh.ucd>"
             << " <merged-mesh.ucd>"
             << endl;
        return 0;
    }
    assert(argc == 4);

    first_filename = argv[1];
    second_filename = argv[2];
    merged_filename = argv[3];

    // read the meshes we want to merge
    cout << "Reading " << first_filename << endl;
    first_mesh.read_ucd(first_filename.c_str());
    cout << "Reading " << second_filename << endl;
    second_mesh.read_ucd(second_filename.c_str());

    // merge them
    cout << "Merging meshes..." << endl;
    union_merge(merged_mesh, first_mesh, second_mesh);

    // write out the result
    cout << "Writing " << merged_filename << endl;
    merged_mesh.write_ucd(merged_filename.c_str());

    return 0;
}
