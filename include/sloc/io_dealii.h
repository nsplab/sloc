#ifndef IO_DEALII_H
#define IO_DEALII_H

#include <string>
#include <vector>
#include <fstream>

#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>

namespace sloc
{
// ----------------------------------------------------------------------------


//
// Output
//

void write_points(const char *filename, const std::vector<dealii::Point<3> >& points);
void write_vector(const char *filename, const dealii::Vector<double>& vec);
void write_matrix(const char *filename, const dealii::FullMatrix<double>& mat);
void write_triangulation(const char *filename, const dealii::Triangulation<2,3>& tria);


//
// Input
//

template <int dim, int spacedim>
void read_ucd_mesh(const char *filename, dealii::Triangulation<dim,spacedim>& tria, bool debug=false)
{
    using namespace std;
    using namespace dealii;

    deallog << "read_ucd_mesh() filename=" << filename << endl;

    GridIn<dim,spacedim> grid_in;
    grid_in.attach_triangulation(tria);

    ifstream infile(filename);
    grid_in.read_ucd(infile);

    if (debug)
    {
        cout << "  Number of active cells: " << tria.n_active_cells() << endl;
        cout << "  Total number of cells: " << tria.n_cells() << endl;
    }
}


// ----------------------------------------------------------------------------
}

#endif // IO_DEALII_H
