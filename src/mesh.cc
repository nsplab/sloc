#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include "mesh.h"
#include "io_ucd.h"

using namespace std;

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

void Mesh::get_point(int n, double *point)
{
    for (int i = 0; i < _ndim; i++)
        point[i] = _pts[_ndim * n + i];
}

void Mesh::get_cell(int e, long *cell)
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


// ----------------------------------------------------------------------------

void ucd_write(const char *filename, Mesh& mesh)
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

void ucd_read(const char *filename, Mesh& mesh)
{
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

// EOF
