#ifndef MESH_H
#define MESH_H

#include <valarray>

class Mesh
{
public:
    Mesh();
    ~Mesh();

    int n_dim() const;
    int n_cells() const;
    int n_points() const;


    void init_points(int npts, int ndim);
    void clear_points();

    void init_cells(int ncells, int ncellnodes);
    void clear_cells();

    void get_point(int n, double *point);
    void get_cell(int e, long *cell);

private:
    int _ndim;
    int _npoints;
    int _ncells;
    int _ncellnodes;
    double *_pts;
    long *_cells;
};

inline int Mesh::n_dim() { return _ndim; }
inline int Mesh::n_points() { return _npoints; }
inline int Mesh::n_cells() { return _ncells; }
inline int Mesh::n_cell_nodes() { return _ncellnodes; }

#endif
