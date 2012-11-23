/*
 * Print mesh information:
 *  - Number of nodes
 *  - Number of elements
 *  - Bounding box
 */

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <sloc/colors.h>

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <sloc/io_stl.h>
#include <getfem/getfem_mesh.h>

using namespace std;
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------

template <typename T>
ostream& operator<<(ostream& os, vcg::Point3<T>& p)
{
    os << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, vcg::Box3<T>& b)
{
    os << "Box3(min=" << b.min << ", max=" << b.max << ")";
    return os;
}

void calc_bbox(sloc::STL_File& stl, vcg::Box3d& bbox)
{
    int i,j;
    float x,y,z;
    for (i = 0; i < stl.n_facets(); ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            stl.get_facet_vertex(i, j, x, y, z);
            bbox.Add(vcg::Point3d(x,y,z));
        }
    }
}

void calc_bbox(getfem::mesh& mesh, vcg::Box3d& bbox)
{
    for (dal::bv_visitor ip(mesh.points_index()); !ip.finished(); ++ip)
    {
        bgeot::base_node p = mesh.points()[ip];
        bbox.Add(vcg::Point3d(p[0], p[1], p[2]));
    }
}


// ----------------------------------------------------------------------------

void print_stl_info(string filename)
{
    cerr << "Calling print_stl_info() with filename=" << filename << endl;

    sloc::STL_File stl;
    stl.read(filename.c_str());

    vcg::Box3d bbox;
    calc_bbox(stl, bbox);

    vcg::Point3d center = bbox.Center();
    vcg::Point3d dim = bbox.Dim();

    cout << ANSI_BLUE << "path   = " << ANSI_CYAN << filename << endl;
    cout << ANSI_BLUE << "cells  = " << ANSI_CYAN << stl.n_facets() << endl;
    cout << ANSI_BLUE << "points = " << ANSI_CYAN << 3 * stl.n_facets() << endl;
    cout << ANSI_BLUE << "bbox   = " << ANSI_CYAN << bbox << endl;
    cout << ANSI_BLUE << "center = " << ANSI_CYAN << center << endl;
    cout << ANSI_BLUE << "dim    = " << ANSI_CYAN << dim << endl;
    cout << ANSI_RESET;
}

void print_getfem_info(string filename)
{
    cerr << "Calling print_getfem_info() with filename=" << filename << endl;

    getfem::mesh mesh;
    mesh.read_from_file(filename);

    vcg::Box3d bbox;
    calc_bbox(mesh, bbox);

    vcg::Point3d center = bbox.Center();
    vcg::Point3d dim = bbox.Dim();

    cout << ANSI_BLUE << "path     = " << ANSI_CYAN << filename << endl;
    cout << ANSI_BLUE << "convexes = " << ANSI_CYAN << mesh.nb_convex() << endl;
    cout << ANSI_BLUE << "points   = " << ANSI_CYAN << mesh.nb_points() << endl;
    cout << ANSI_BLUE << "bbox     = " << ANSI_CYAN << bbox << endl;
    cout << ANSI_BLUE << "center   = " << ANSI_CYAN << center << endl;
    cout << ANSI_BLUE << "dim      = " << ANSI_CYAN << dim << endl;
    cout << ANSI_RESET;
}

int print_info(fs::path meshpath)
{
    string ext = meshpath.extension().string();

    if (ext == ".mesh")
    {
        print_getfem_info(meshpath.string());
    }
    else if (ext == ".stl")
    {
        print_stl_info(meshpath.string());
    }
    else
    {
        cerr << "File extension " << ANSI_RED << ext << ANSI_RESET << " is not supported!" << endl;
        exit(1);
    }

    return 0;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " MESHFILE" << endl;
        return 0;
    }

    fs::path meshpath(argv[1]);
    if (!fs::exists(meshpath))
    {
        cerr << "File " << ANSI_RED << meshpath.string() << ANSI_RESET << " does not exist!" << endl;
        return 1;
    }

    if (!meshpath.has_extension())
    {
        cerr << "File " << ANSI_RED << meshpath.string() << ANSI_RESET << " has no file extension!" << endl;
        return 1;
    }

    print_info(meshpath);

    return 0;
}
