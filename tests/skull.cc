/* skull.cc
 *
 * Program to convert Bone2.stl and Artery2.stl into .ucd files
 * In the process, we also need to convert the corresponding
 * triangle mesh into a quadrilateral mesh.
 */
#include <iostream>
#include <ctime>
#include "mesh.h"
#include "tri2quad.h"

using namespace std;
using namespace sloc;

int main(void)
{
    Mesh tmesh, qmesh;

    // process Bone2.stl
    tmesh.read_stl("tmp/Bone2.stl");
    tri2quad(tmesh, qmesh);
    tmesh.write_ucd("tmp/Bone2.tri.ucd");
    qmesh.write_ucd("tmp/Bone2.quad.ucd");

    // process Artery2.stl
    tmesh.read_stl("tmp/Artery2.stl");
    tri2quad(tmesh, qmesh);
    tmesh.write_ucd("tmp/Artery2.tri.ucd");
    qmesh.write_ucd("tmp/Artery2.quad.ucd");

    return 0;
}
