#ifndef SLOC_MESH_H
#define SLOC_MESH_H

#include <valarray>

namespace sloc
{
// ----------------------------------------------------------------------------

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

    void get_point(int n, double *point) const;
    void get_cell(int e, long *cell) const;
    void get_mat(int e, int& mat) const;

    void set_point(int n, double *point);
    void set_cell(int e, long *cell);
    void set_mat(int e, int mat);

    int write(const char *filename);
    int read(const char *filename);

    void write_ucd(const char *filename);
    void read_ucd(const char *filename);
    void read_stl(const char *filename);

private:
    int _ndim;
    int _npoints;
    int _ncells;
    int _ncellnodes;
    double *_pts;
    long *_cells;
    int *_mat;
};

inline int Mesh::n_dim() const { return _ndim; }
inline int Mesh::n_points() const { return _npoints; }
inline int Mesh::n_cells() const { return _ncells; }
inline int Mesh::n_cell_nodes() const { return _ncellnodes; }

void ucd_write(Mesh& mesh, const char *filename);
void stl_write(Mesh& mesh, const char *filename);

void ucd_read(Mesh& mesh, const char *filename);
void stl_read(Mesh& mesh, const char *filename);


// ----------------------------------------------------------------------------
}
#endif
