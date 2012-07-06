/*
 * Test whether given mesh can be loaded with Deal.II
 */

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include "color_utils.h"
//#include "mesh.h"
#include "io_dealii.h"

namespace fs = boost::filesystem;
using namespace std;

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        cout << "Usage: "
             << ANSI_RED << argv[0] << " MESHFILE"
             << ANSI_RESET << endl;
        return 0;
    }

    fs::path mesh_path(argv[1]);
    if (!fs::exists(mesh_path))
    {
        cerr << "File "
             << ANSI_RED << mesh_path.filename()
             << ANSI_RESET << " does not exist!" << endl;
        return 1;
    }

    const bool debug = true;
    dealii::Triangulation<2,3> tria;
    sloc::read_ucd_mesh(mesh_path.c_str(), tria, debug);

    return 0;
}
