/* stl_to_ucd.cc
 *
 * Convert stl file to a ucd mesh.
 *
 */

#include <iostream>
#include <string>
#include <cassert>
#include "mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
    sloc::Mesh mesh;
    string infile, outfile;

    if (argc == 1)
    {
        cout << "Usage: " << argv[0]
             << " <input-mesh.stl> <output-mesh.ucd>"
             << endl;
        return 0;
    }
    assert(argc == 3);

    infile = argv[1];
    outfile = argv[2];

    cout << "Reading " << infile << endl;
    mesh.read_stl(infile.c_str());

    cout << "Writing " << outfile << endl;
    mesh.write_ucd(outfile.c_str());

    return 0;
}
