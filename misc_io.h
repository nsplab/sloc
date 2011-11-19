#ifndef MISC_IO_H
#define MISC_IO_H

#include <string>
#include <vector>
#include <deal.II/grid/tria.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>


// for reading ucd meshes
template <int dim, int spacedim>
void read_ucd_mesh(const std::string &filename, dealii::Triangulation<dim,spacedim> &tria);

// for reading cubit's exodus meshes
void read_exo_surface_mesh(const std::string &filename, dealii::Triangulation<2,3> &tria);

// useful for debugging
void write_triangulation(const std::string &filename, const dealii::Triangulation<2,3> &tria);
void write_points(std::ofstream &out, const std::vector<dealii::Point<3> > &points);
void write_points(const std::string &filename, const std::vector<dealii::Point<3> > &points);
void write_vector(const std::string &filename, const dealii::Vector<double> &vec);
void write_matrix(const std::string &filename, const dealii::FullMatrix<double> &mat);

#endif // MISC_IO_H
