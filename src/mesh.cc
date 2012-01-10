#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cstring>
#include "mesh.h"
#include "io_ucd.h"
#include "io_stl.h"
#include "point_cloud.h"
#include "progress_timer.h"

using namespace std;
using namespace sloc;

Mesh::Mesh()
{
    _ndim = 0;
    _npoints = 0;
    _ncells = 0;
    _ncellnodes = 0;
    _pts = 0;
    _cells = 0;
}

Mesh::~Mesh()
{
    clear_points();
    clear_cells();
}

void Mesh::init_points(int npts, int ndim)
{
    clear_points();
    _npoints = npts;
    _ndim = ndim;
    _pts = new double[_npoints * _ndim];
}

void Mesh::init_cells(int ncells, int ncellnodes)
{
    clear_cells();
    _ncells = ncells;
    _ncellnodes = ncellnodes;
    _cells = new long[_ncells * _ncellnodes];
}

void Mesh::clear_points()
{
    if (_pts != 0) delete [] _pts;
    _ndim = 0;
    _npoints = 0;
    _pts = 0;
}

void Mesh::clear_cells()
{
    if (_cells != 0) delete [] _cells;
    _ncells = 0;
    _ncellnodes = 0;
    _cells = 0;
}

void Mesh::get_point(int n, double *point) const
{
    for (int i = 0; i < _ndim; i++)
        point[i] = _pts[_ndim * n + i];
}

void Mesh::get_cell(int e, long *cell) const
{
    for (int i = 0; i < _ncellnodes; i++)
        cell[i] = _cells[_ncellnodes * e + i];
}

void Mesh::set_point(int n, double *point)
{
    for (int i = 0; i < _ndim; i++)
        _pts[_ndim * n + i] = point[i];
}

void Mesh::set_cell(int e, long *cell)
{
    for (int i = 0; i < _ncellnodes; i++)
        _cells[_ncellnodes * e + i] = cell[i];
}

void Mesh::write_ucd(const char *filename)
{
    cout << "Writing to " << filename << endl;
    ucd_write(filename, *this);
}

void Mesh::read_ucd(const char *filename)
{
    cout << "Reading from " << filename << endl;
    ucd_read(filename, *this);
}

void Mesh::read_stl(const char *filename)
{
    cout << "Reading from " << filename << endl;
    stl_read(filename, *this);
}

// ----------------------------------------------------------------------------

void sloc::ucd_write(const char *filename, Mesh& mesh)
{
    int e,n,i;
    ofstream file;
    file.open(filename);

    file << "# ucd file" << endl;

    file << mesh.n_points() << " "
         << mesh.n_cells() << " "
         << "0 0 0"
         << endl;

    // Write out nodes
    for (n = 0; n < mesh.n_points(); n++)
    {
        double point[3];
        mesh.get_point(n, point);
        file << n
             << " " << point[0]
             << " " << point[1]
             << " " << point[2]
             << endl;
    }

    // TODO: generalize for other cell types
    char celltype[10];
    if (mesh.n_cell_nodes() == 3)
        strcpy(celltype, "tri");
    else if (mesh.n_cell_nodes() == 4)
        strcpy(celltype, "quad");
    else
        cout << "bad cell type!\n";

    // Write out elements
    long *cell = new long[mesh.n_cell_nodes()];
    for (e = 0; e < mesh.n_cells(); e++)
    {
        mesh.get_cell(e, cell);

        // TODO: get the material id's from mesh object.
        // for now, just assume it's always 1.
        long material_id = 1;

        file << e
             << " " << material_id
             << " " << celltype;

        for (i = 0; i < mesh.n_cell_nodes(); i++)
        {
            file << " " << cell[i];
        }
        file << endl;
    }
    delete [] cell;

    file.close();
}

void sloc::ucd_read(const char *filename, Mesh& mesh)
{
    cout << "sloc::ucd_read()\n";
    cout << "  filename = " << filename << endl;
    UCD_File ucd;
    ucd.read(filename);

    // copy nodes
    mesh.init_points(ucd.num_nodes, 3);
    for (int n = 0; n < mesh.n_points(); n++)
    {
        double point[3];
        point[0] = ucd._nodes[3*n+0];
        point[1] = ucd._nodes[3*n+1];
        point[2] = ucd._nodes[3*n+2];
        mesh.set_point(n, point);
    }

    // copy cells
    mesh.init_cells(ucd.num_cells, ucd._cells[0]->num_vertices());
    long *cell = new long[mesh.n_cell_nodes()];
    for (int e = 0; e < mesh.n_cells(); e++)
    {
        UCD_Cell *c = ucd._cells[e];
        for (int i = 0; i < mesh.n_cell_nodes(); i++)
            cell[i] = c->cell_verts[i];
        mesh.set_cell(e, cell);
    }
    delete [] cell;

    ucd.clear();
}

void sloc::stl_read(const char *filename, Mesh& tmesh)
{
    //
    // Given the data in the stl file, we need to create
    // a mesh object with all duplicate points removed.
    //

    cout << "sloc::stl_read()\n";
    cout << "  filename = " << filename << endl;

    STL_File stl;
    stl.read(filename);

    // duplicate point elimination takes a while, so we need
    // something to print out periodic progress reports to stdout.
    time_t t0, t1;
    ProgressTimer timer;

    // loop indices
    int e, n;

    // our triangular mesh gets its cells from the STL facets.
    int ncells = stl.n_facets();
    int npts = ncells * 3;

    cout << "  stl points = " << npts << endl;
    cout << "  stl facets = " << ncells << endl;

    long *ids = new long[npts];
    double *pts = new double[npts * 3];

    //
    // to determine the list of unique points, we load them
    // into a point cloud structure that can detect duplicates
    // within a given tolerance.
    //
    PointCloud points;
    points.set_tolerance(1e-8);

    cout << "  Building point cloud...\n";
    time(&t0);

    cout << timer.header("facets");
    timer.start(ncells);

    for (e = 0; e < ncells; e++)
    {
        points.add(stl.va[3*e+0], stl.va[3*e+1], stl.va[3*e+2], &ids[3*e+0]);
        points.add(stl.vb[3*e+0], stl.vb[3*e+1], stl.vb[3*e+2], &ids[3*e+1]);
        points.add(stl.vc[3*e+0], stl.vc[3*e+1], stl.vc[3*e+2], &ids[3*e+2]);
        if (e % 1000 == 0)
            cout << timer.update(e);
    }
    cout << timer.update(ncells) << endl;

    time(&t1);
    cout << "  Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";

    // now, we can query the point cloud for its count,
    // and use that to initialize the mesh
    tmesh.init_points(points.n_points(), 3);
    tmesh.init_cells(ncells, 3);

    cout << "  mesh points = " << tmesh.n_points() << endl;
    cout << "  mesh cells  = " << tmesh.n_cells() << endl;

    // load points from point cloud into mesh
    for (n = 0; n < tmesh.n_points(); n++)
    {
        double pt[3];
        points.get_point(n, pt);
        tmesh.set_point(n, pt);
    }

    // load cells into mesh
    for (e = 0; e < tmesh.n_cells(); e++)
    {
        tmesh.set_cell(e, &ids[3*e]);
    }

    delete [] pts;
    delete [] ids;
}

// EOF
