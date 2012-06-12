/* tri_to_quad.cc
 *
 * Convert triangular mesh to a quadrilateral mesh.
 *
 */

#include <iostream>
#include <string>
#include <cassert>
#include "mesh.h"
#include "geometry_utils.h"

using namespace std;

int main(int argc, char *argv[])
{
    sloc::Mesh tmesh, qmesh;
    string infile, outfile;

    if (argc == 1)
    {
        cout << "Usage: " << argv[0]
             << " <input-tri-mesh.ucd> <output-quad-mesh.ucd>"
             << endl;
        return 0;
    }
    assert(argc == 3);

    infile = argv[1];
    outfile = argv[2];

    cout << "Reading " << infile << endl;
    tmesh.read_ucd(infile.c_str());

    cout << "Converting triangles to quads..." << endl;
    sloc::tri2quad(tmesh, qmesh, true);

    cout << "Writing " << outfile << endl;
    qmesh.write_ucd(outfile.c_str());

    return 0;
}
