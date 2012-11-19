/*
 * Save GetFEM native meshes to VTK format for visualization.
 *
 */
#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_export.h>

using namespace std;
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << argv[0] << " input.mesh output.vtk" << endl;
        return 0;
    }

    string infile = argv[1];
    string outfile = argv[2];

    if (!fs::exists(infile))
    {
        cerr << "File does not exist: " << infile << endl;
        return 1;
    }

    getfem::mesh m;
    m.read_from_file(infile);
    cerr << "Read " << infile << endl;

    getfem::vtk_export exp(outfile.c_str(), true);
    exp.exporting(m);
    exp.write_mesh();
    cerr << "Wrote " << outfile << endl;

    return 0;
}
