/*
 * mesh_info.cc
 *
 * Print out mesh information
 */

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include "mesh.h"
#include "color_utils.h"

namespace fs = boost::filesystem;
using namespace std;

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " MESHFILE" << endl;
        return 0;
    }

    fs::path mesh_path(argv[1]);

    if (!fs::exists(mesh_path))
    {
        cerr << "File " << mesh_path.filename() << " does not exist!" << endl;
        return 1;
    }

    if (!mesh_path.has_extension())
    {
        cerr << "Path " << mesh_path << " has no file extension!" << endl;
        return 2;
    }

    string extension = mesh_path.extension().string();

    sloc::Mesh mesh;
    if (extension == ".ucd")
        mesh.read_ucd(mesh_path.c_str());
    else if (extension == ".stl")
        mesh.read_stl(mesh_path.c_str());
    else
    {
        cerr << "Unsupported file extension '" << extension << "'" << endl;
        return 3;
    }

    cout << ANSI_BLUE << "path   = " << ANSI_CYAN << mesh_path << endl;
    cout << ANSI_BLUE << "points = " << ANSI_CYAN << mesh.n_points() << " x " << mesh.n_dim() << endl;
    cout << ANSI_BLUE << "cells  = " << ANSI_CYAN << mesh.n_cells() << " x " << mesh.n_cell_nodes() << endl;
    cout << ANSI_RESET;

    return 0;
}

