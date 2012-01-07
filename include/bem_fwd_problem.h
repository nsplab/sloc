#ifndef SLOC_BEM_FWD_H
#define SLOC_BEM_FWD_H

/* deal.II includes */

// for handling triangulations
#include <deal.II/grid/tria.h>

// linear algebra includes
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>

// for managing the distribution and number of the degrees of freedom
// (i.e., associating dofs with vertices, lines, and cells).
#include <deal.II/dofs/dof_handler.h>

// for the Lagrange finite elements we'll use
#include <deal.II/fe/fe_q.h>

// for the Gauss-Legendre quadrature rules we'll be using
#include <deal.II/base/quadrature_lib.h>

// for mapping reference element to arbitrary cell
#include <deal.II/fe/mapping_q.h>

// various utility functions
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

/* other includes */
#include "dipole_sources.h"
#include "material_data.h"

namespace sloc
{
// ----------------------------------------------------------------------------

class BEM_ForwardProblem
{
public:

    BEM_ForwardProblem();
    ~BEM_ForwardProblem();

    void run();

private:

    void configure();
    void assemble_system();
    void solve_system();
    void output_results();

    void compute_general_solution();
    void compute_area();

private:

    // XXX: what is this sigma0 everyone keeps talking about
    static const double sigma0 = 1.0;

    // location and components of current dipole source
    DipoleSources dipole_sources;
    
    // conductivity values across interface surfaces (sigma_int, sigma_ext)
    MaterialData material_data;

    // triangulation of domain: use 2d cells embedded in 3d space
    dealii::Triangulation<2,3> tria;

    // use a 2d Lagrange finite element embedded in 3d space
    dealii::FE_Q<2,3> fe;

    // degrees of freedom for 2d elements embedded in 3d space
    dealii::DoFHandler<2,3> dh;

    // the map function between 2d reference element and its cell in 3d space
    dealii::MappingQ<2,3> mapping;

    // the 2d quadrature rule to use when integrating our 2d reference element
    dealii::Quadrature<2> quadrature;

    // the system matrix (I-C), and the rhs vector g.
    dealii::FullMatrix<double> system_matrix;
    dealii::Vector<double> system_rhs;

    // the solution $\phi$
    dealii::Vector<double> phi;

    // parameters for controlling linear solver (tolerance, logging, etc.)
    dealii::SolverControl solver_control;

    // for measuring wallclock time
    dealii::Timer timer;
};

// ----------------------------------------------------------------------------
}
#endif
