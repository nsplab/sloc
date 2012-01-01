// for handling triangulations
#include <deal.II/grid/tria.h>

// for looping over cells and faces
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

// input of grids in various formats
#include <deal.II/grid/grid_in.h>

// output of grids in various graphics formats
#include <deal.II/grid/grid_out.h>

// standard includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "io_dealii.h"

// namespaces
using namespace std;
using namespace dealii;

// ----------------------------------------------------------------------------

int main()
{
    deallog << "main()" << endl;

    Triangulation<3> tria;
    sloc::read_ucd_mesh("doublesphere.ucd", tria);

    Triangulation<3>::active_cell_iterator
      cell,
      endc = tria.end();

    int n = GeometryInfo<3>::vertices_per_cell;

    ofstream out;
    out.open("doublesphere.mat");

    for (cell = tria.begin_active(); cell != endc; ++cell)
    {
      Point<3> center(0,0,0);

      for(int v = 0; v < n; ++v)
      {
        center += cell->vertex(v);
      }
      center /= n;

      double distance2 = center.square();
      if (distance2 <= 4)
      {
        out << "1" << endl;
      }
      else
      {
        out << "2" << endl;
      }
    }
    out.close();
    cout << "Wrote material id to file doublesphere.mat" << endl;
    return 0;
}
