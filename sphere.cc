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

#include "misc_io.h"

// namespaces
using namespace std;
using namespace dealii;

// ----------------------------------------------------------------------------

int main()
{
    deallog << "main()" << endl;

    // surface triangulation
    Triangulation<2,3> tria2;
    read_ucd_mesh("tmp/sphere.surf.ucd", tria2);

    // volume triangulation
    Triangulation<3> tria3;
    read_ucd_mesh("tmp/sphere.ucd", tria3);

    return 0;
}
