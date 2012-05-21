/* bone_artery_stl_to_ucd.cc
 *
 * Program to convert Bone2.stl and Artery2.stl into .ucd files
 * In the process, we also need to convert the corresponding
 * triangle mesh into a quadrilateral mesh.
 *
 * Use this program to convert:
 *   Artery2.stl -> Artery2.quad.ucd
 *   Bone2.stl -> Bone2.quad.ucd
 */
#include <iostream>
#include <ctime>
#include "mesh.h"
#include "geometry_utils.h"
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace sloc;

int main(int argc, char *argv[])
{
    Mesh tmesh, qmesh;

    if (argc == 1)
    {
        cout << "Usage: " << argv[0] << " file.stl" << endl;
        return 0;
    }

    string infile, prefix;
    string outfile1, outfile2;

    if (boost::algorithm::ends_with(argv[1], ".stl"))
    {
        infile = argv[1];
        prefix = infile.substr(0, infile.size()-4);
    }
    else
    {
        prefix = argv[1];
        infile = prefix + ".stl";
    }

    outfile1 = prefix + ".tri.ucd";
    outfile2 = prefix + ".quad.ucd";

    cout << "Reading " << infile << endl;
    tmesh.read_stl(infile.c_str());

    cout << "Converting triangles to quads..." << endl;
    tri2quad(tmesh, qmesh);

    cout << "Writing " << outfile1 << endl;
    tmesh.write_ucd(outfile1.c_str());

    cout << "Writing " << outfile2 << endl;
    qmesh.write_ucd(outfile2.c_str());

    return 0;
}
