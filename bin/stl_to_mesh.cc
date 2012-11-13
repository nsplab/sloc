/*
 * Convert STL file into format suited for GetFEM
 */

#include <iostream>
#include <string>
#include <sloc/io_stl.h>
#include <sloc/io_getfem.h>
#include <getfem/getfem_mesh.h>

using namespace std;
int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cout << "Usage: " << argv[0] << " <infile.stl> <outfile.mesh>" << endl;
        return 0;
    }

    string infile = argv[1];
    string outfile = argv[2];

    getfem::mesh mesh;
    sloc::read_stl(mesh, infile.c_str());
    mesh.write_to_file(outfile);

    return 0;
}
