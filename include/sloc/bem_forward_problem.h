#ifndef BEM_FORWARD_PROBLEM_H
#define BEM_FORWARD_PROBLEM_H


// mesh, fem, and quadrature methods
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_mesh_im.h>
//#include <getfem/getfem_import.h>
#include <getfem/bgeot_geometric_trans.h>
//#include <getfem/getfem_integration.h>


// linear algebra includes
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>

// various utility classes and functions
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

// other includes
#include <fstream>
#include <sloc/dipole_sources.h>
#include <sloc/material_data.h>


namespace sloc {
// ----------------------------------------------------------------------------

class BEM_Forward_Problem
{
public:

    class Parameters
    {
    public:
        static void declare_parameters(dealii::ParameterHandler& prm);
        void get_parameters(dealii::ParameterHandler& prm);
        void print() const;
    public:
        std::string surface_mesh;
        std::string surface_mesh_materials;
        std::string material_data;
        std::string dipole_sources;
        std::string output_vtk;
        std::string output_phi;
        std::string logfile;
        bool debug;
        bool verbose;
    };

public:
    BEM_Forward_Problem(const Parameters& parameters);
    ~BEM_Forward_Problem();

    void run();

protected:

    void configure();
    void assemble_system();
    void solve_system();
    void output_results();

    void compute_general_solution();
    void compute_area();

protected:

    // this "parameters" object encapsulates our access to any parameter
    // we might need during the configuration step
    Parameters parameters;

    // location and components of the current dipole sources
    DipoleSources dipole_sources;

    // conductivity values across interface surfaces (sigma_int and sigma_ext)
    MaterialData material_data;

    // finite element method on mesh: FEM_PK(2,1)
    getfem::pfem pf;

    // the map function between the 2d reference element and its cell in 3d space
    bgeot::pgeometric_trans pgt;

    // integration method on mesh
    getfem::pintegration_method pim;

    // surface mesh -- 2d cells embedded in 3d space
    getfem::mesh surface_mesh;
    getfem::mesh_fem mf;
    getfem::mesh_im mim;

    // the system matrix (I-C), and the rhs vector g
    dealii::FullMatrix<double> system_matrix;
    dealii::Vector<double> system_rhs;

    // the solution $\phi$
    dealii::Vector<double> phi;

    // parameters for controlling linear solver (tolerance, logging, etc.)
    dealii::SolverControl solver_control;

    // for measuring wallclock time
    dealii::Timer timer;

    // output behavior flags
    bool verbose;
    bool debug;

    // log stream for debug mode
    std::ofstream log;

};

// ----------------------------------------------------------------------------
}

#endif
