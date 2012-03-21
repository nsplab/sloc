/* head.cc
 *
 * Merge Bone2.quad.ucd and Artery2_1.quad.ucd into a single file called head.ucd
 *
 */

#include <string>
#include <iostream>
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

int main(void)
{
    string bone_file = "tmp/Bone2_1.quad.ucd";
    string artery_file = "tmp/Artery2_1.quad.ucd";
    string merged_file = "tmp/head.ucd";

    sloc::Mesh merged_mesh, bone_mesh, artery_mesh;

    // read the meshes we want to merge
    bone_mesh.read_ucd(bone_file.c_str());
    artery_mesh.read_ucd(artery_file.c_str());

    // merge them
    union_merge(merged_mesh, bone_mesh, artery_mesh);

    // write out the result
    merged_mesh.write_ucd(merged_file.c_str());

    return 0;
}
