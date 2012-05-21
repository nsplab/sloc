/* set_material_id.cc - set the material ID of a UCD mesh
 *
 */

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include "mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        cout << "Usage: " << argv[0]
             << " <new-mat-id> <in-mesh.ucd> <out-mesh.ucd>"
             << endl;
        return 0;
    }
    assert(argc == 4);

    int mat_id = boost::lexical_cast<int>(argv[1]);
    string infile = argv[2];
    string outfile = argv[3];

    sloc::Mesh mesh;

    cout << "Reading " << infile << endl;
    mesh.read_ucd(infile.c_str());

    cout << "Setting material id of every cell to " << mat_id << endl;
    for (int e = 0; e < mesh.n_cells(); e++)
        mesh.set_mat(e, mat_id);

    cout << "Writing " << outfile << endl;
    mesh.write_ucd(outfile.c_str());

    return 0;
}
