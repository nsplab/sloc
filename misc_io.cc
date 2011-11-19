#include "misc_io.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include <netcdfcpp.h>

using namespace std;
using namespace dealii;

template <int dim, int spacedim>
void read_ucd_mesh(const std::string &filename, Triangulation<dim,spacedim> &tria)
{
    deallog << "read_ucd_mesh()" << endl;

    GridIn<dim,spacedim> grid_in;
    grid_in.attach_triangulation(tria);

    std::ifstream input_file(filename.c_str());
    grid_in.read_ucd(input_file);

    if (true)
    {
        cout << "Number of active cells: " << tria.n_active_cells() << endl;
        cout << "Total number of cells : " << tria.n_cells() << endl;
    }
}

// Here's a function to read mesh from a NetCDF file
// For now, we just assume our information is stored in block 1
void read_exo_surface_mesh(const std::string &filename, Triangulation<2,3> &tria)
{
    deallog << "read_exo_surface_mesh()" << endl;

    const bool output = false;

    NcFile nc(filename.c_str());


    NcDim *elements_dim = nc.get_dim("num_el_in_blk1");
    AssertThrow(elements_dim->is_valid(), ExcIO());
    const unsigned int n_cells = elements_dim->size();
    if (output)
        cout << "n_cells=" << n_cells << endl;


    // num_nod_per_el1
    NcDim *num_nod_per_el1 = nc.get_dim("num_nod_per_el1");
    AssertThrow(num_nod_per_el1->is_valid(), ExcIO());
    const unsigned int vertices_per_cell = num_nod_per_el1->size();
    AssertThrow(vertices_per_cell==GeometryInfo<2>::vertices_per_cell, ExcIO());


    //
    NcVar *vertex_indices_var = nc.get_var("connect1");
    AssertThrow(vertex_indices_var->is_valid(), ExcIO());
    AssertThrow(vertex_indices_var->num_dims()==2, ExcIO());
    AssertThrow(static_cast<unsigned int>(vertex_indices_var->get_dim(0)->size())==n_cells, ExcIO());
    AssertThrow(static_cast<unsigned int>(vertex_indices_var->get_dim(1)->size())==vertices_per_cell, ExcIO());

    vector<int> vertex_indices(n_cells * vertices_per_cell);
    vertex_indices_var->get(&*vertex_indices.begin(), n_cells, vertices_per_cell);

    if (output)
    {
        cout << "vertex_indices:" << endl;
        for (unsigned int i=0, v=0; i < n_cells; ++i)
        {
            for (unsigned int j=0; j < vertices_per_cell; ++j)
                cout << vertex_indices[v++] << " ";
            cout << endl;
        }
    }


    //
    NcDim *vertices_dim = nc.get_dim("num_nodes");
    AssertThrow(vertices_dim->is_valid(), ExcIO());
    const unsigned int n_vertices = vertices_dim->size();
    if (output)
        cout << "n_vertices=" << n_vertices << endl;

    NcVar *coordx = nc.get_var("coordx");
    NcVar *coordy = nc.get_var("coordy");
    NcVar *coordz = nc.get_var("coordz");
    AssertThrow(coordx->is_valid(), ExcIO());
    AssertThrow(coordy->is_valid(), ExcIO());
    AssertThrow(coordz->is_valid(), ExcIO());
    AssertThrow(coordx->num_dims()==1, ExcIO());
    AssertThrow(coordy->num_dims()==1, ExcIO());
    AssertThrow(coordz->num_dims()==1, ExcIO());

    vector<vector<double> > point_values(3, vector<double>(n_vertices));
    coordx->get(&*point_values[0].begin(), n_vertices);
    coordy->get(&*point_values[1].begin(), n_vertices);
    coordz->get(&*point_values[2].begin(), n_vertices);

    vector<Point<3> > vertices(n_vertices);
    for (unsigned int i = 0; i < n_vertices; i++)
    {
        vertices[i](0) = point_values[0][i];
        vertices[i](1) = point_values[1][i];
        vertices[i](2) = point_values[2][i];
    }

    //
    vector<CellData<2> > cells(n_cells);
    for (unsigned int cell=0; cell < n_cells; ++cell)
    {
        for (unsigned int i=0; i < vertices_per_cell; ++i)
            cells[cell].vertices[i] = vertex_indices[cell*vertices_per_cell + i] - 1;
    }

    SubCellData subcelldata;
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    GridReordering<2,3>::reorder_cells(cells);
    tria.create_triangulation_compatibility(vertices, cells, subcelldata);
}

void write_triangulation(const std::string &filename, const dealii::Triangulation<2,3> &tria)
{
    deallog << "write_triangulation()" << endl;
    ofstream out;
    out.open(filename.c_str());
    GridOut grid_out;
    grid_out.write_ucd(tria, out);
    out.close();
}

void write_points(std::ofstream &out, const std::vector<dealii::Point<3> > &points)
{
    for (unsigned int i = 0; i < points.size(); ++i)
        out << points[i](0) << " "
            << points[i](1) << " "
            << points[i](2) << endl;
}

void write_points(const std::string &filename, const std::vector<dealii::Point<3> > &points)
{
    deallog << "write_points()" << endl;
    ofstream out;
    out.open(filename.c_str());
    write_points(out, points);
    out.close();
}


void write_vector(const std::string &filename, const dealii::Vector<double> &vec)
{
    deallog << "write_vector()" << endl;
    ofstream out;
    out.open(filename.c_str());
    for (unsigned int i = 0; i < vec.size(); ++i)
        out << vec(i) << endl;
    out.close();
}

void write_matrix(const std::string &filename, const dealii::FullMatrix<double> &mat)
{
    deallog << "write_matrix()" << endl;
    ofstream out;
    out.open(filename.c_str());
    for (unsigned int i = 0; i < mat.m(); ++i)
    {
        for (unsigned int j = 0; j < mat.n(); ++j)
            out << mat(i,j) << " ";
        out << endl;
    }
    out.close();
}
