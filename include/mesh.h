#ifndef MESH_H
#define MESH_H

#include <valarray>

namespace sloc
{

class Mesh
{
public:
    Mesh();
    ~Mesh();

    int n_dim() const;
    int n_points() const;
    int n_cells() const;
    int n_cell_nodes() const;

    void init_points(int npts, int ndim);
    void clear_points();

    void init_cells(int ncells, int ncellnodes);
    void clear_cells();

    void get_point(int n, double *point);
    void get_cell(int e, long *cell);

    void set_point(int n, double *point);
    void set_cell(int e, long *cell);

    void write_ucd(const char *filename);


private:
    int _ndim;
    int _npoints;
    int _ncells;
    int _ncellnodes;
    double *_pts;
    long *_cells;
};

inline int Mesh::n_dim() const { return _ndim; }
inline int Mesh::n_points() const { return _npoints; }
inline int Mesh::n_cells() const { return _ncells; }
inline int Mesh::n_cell_nodes() const { return _ncellnodes; }

void ucd_write(const char *filename, Mesh& mesh);
void ucd_read(const char *filename, Mesh& mesh);

}
#endif
