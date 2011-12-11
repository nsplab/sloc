#include "mesh.h"

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
    _npoints = npts;
    _ndim = ndim;
    _pts = new double[_npoints * _ndim];
}

void Mesh::init_cells(int ncells, int ncellnodes)
{
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

