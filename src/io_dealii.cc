#include "sloc/io_dealii.h"

#include <iostream>
#include <deal.II/grid/grid_out.h>

using namespace std;
using namespace dealii;

void sloc::write_points(const char *filename, const vector<dealii::Point<3> >& points)
{
    deallog << "write_points() filename=" << filename << endl;
    ofstream out;
    out.open(filename);
    for (unsigned int i = 0; i < points.size(); ++i)
    {
        const Point<3>& pt = points[i];
        out << pt(0) << " "
            << pt(1) << " "
            << pt(2) << endl;
    }
    out.close();
}

void sloc::write_vector(const char *filename, const dealii::Vector<double>& vec)
{
    deallog << "write_vector() filename=" << filename << endl;
    ofstream out;
    out.open(filename);
    for (unsigned int i = 0; i < vec.size(); ++i)
        out << vec(i) << endl;
    out.close();
}

void sloc::write_matrix(const char *filename, const dealii::FullMatrix<double>& mat)
{
    deallog << "write_matrix() filename=" << filename << endl;
    ofstream out;
    out.open(filename);
    for (unsigned int i = 0; i < mat.m(); ++i)
    {
        for (unsigned int j = 0; j < mat.n(); ++j)
            out << " " << mat(i,j);
        out << endl;
    }
    out.close();
}

void sloc::write_triangulation(const char *filename, const dealii::Triangulation<2,3>& tria)
{
    deallog << "write_triangulation() filename=" << filename << endl;
    ofstream out;
    out.open(filename);
    GridOut grid_out;
    grid_out.write_ucd(tria, out);
    out.close();
}

