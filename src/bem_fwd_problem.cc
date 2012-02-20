#include "bem_fwd_problem.h"

/* additional deal.II includes */

// for looping over cells and faces
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

// linear algebra
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

// provide information about the degrees of freedom local to a cell
#include <deal.II/dofs/dof_accessor.h>

// for the function DoFTools::map_dofs_to_support_points
#include <deal.II/dofs/dof_tools.h>

// for evaluating the finite element shape functions at quadrature points of cell
#include <deal.II/fe/fe_values.h>

// input of grids in various formats
#include <deal.II/grid/grid_in.h>

// output of grids in various graphic formats
#include <deal.II/grid/grid_out.h>

// for writing out the output file given the shape function coefficients from the solution
#include <deal.II/numerics/data_out.h>

// for the function VectorTools::integrate_difference
#include <deal.II/numerics/vectors.h>

// various utilities
#include <deal.II/base/utilities.h>
#include <deal.II/base/table.h>

/* standard includes */
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

/* other includes */
#include "io_dealii.h"
#include "progress_timer.h"

using namespace dealii;
using namespace sloc;

// ----------------------------------------------------------------------------

void BEM_ForwardProblem::Parameters::declare_parameters(ParameterHandler& prm)
{
    prm.declare_entry("dipole_sources", "", Patterns::Anything(), "Filename for current dipole sources data");
    prm.declare_entry("material_data", "", Patterns::Anything(), "Filename for material data");
    prm.declare_entry("surface_mesh", "", Patterns::Anything(), "Filename for surface mesh");
    prm.declare_entry("volume_mesh", "", Patterns::Anything(), "Filename for volume mesh");
    prm.declare_entry("surface_phi", "phi", Patterns::Anything(), "Filename prefix for vtk output file (surface solution)");
    prm.declare_entry("volume_phi", "phi_vol", Patterns::Anything(), "Filename prefix for vtk output file (volume solution)");
    prm.declare_entry("debug_logfile", "debug.log", Patterns::Anything(), "Filename for debugging logfile");
    prm.declare_entry("debug", "false", Patterns::Bool(), "Output debug information");
    prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level");
}

void BEM_ForwardProblem::Parameters::get_parameters(ParameterHandler& prm)
{
    dipole_sources = prm.get("dipole_sources");
    material_data = prm.get("material_data");
    surface_mesh = prm.get("surface_mesh");
    volume_mesh = prm.get("volume_mesh");
    surface_phi = prm.get("surface_phi");
    volume_phi = prm.get("volume_phi");
    debug_logfile = prm.get("debug_logfile");
    debug = prm.get_bool("debug");
    verbose = prm.get_bool("verbose");
}

// ----------------------------------------------------------------------------

BEM_ForwardProblem::BEM_ForwardProblem(const Parameters& parameters)
    : parameters(parameters),   // parameters class
      fe(1),                    // fe degree 1
      dh(tria),                 // attach triangulation to our dof_handler
      mapping(1, true)          // mapping degree 1 (also, use same mapping on interior cells)
{
    timer.start();
}

BEM_ForwardProblem::~BEM_ForwardProblem()
{
}

void BEM_ForwardProblem::run()
{
    deallog << "BEM_ForwardProblem::run() T=" << timer.wall_time() << std::endl;
    configure();
    compute_area();
    assemble_system();
    solve_system();
    output_results();
    //compute_general_solution();
    deallog << "last timer = " << timer.wall_time() << std::endl;
}

// ----------------------------------------------------------------------------

void BEM_ForwardProblem::configure()
{
    deallog << "BEM_ForwardProblem::configure() T=" << timer.wall_time() << std::endl;

    // configure the dipole sources
    Assert(!parameters.dipole_sources.empty(), ExcEmptyObject());
    dipole_sources.read(parameters.dipole_sources.c_str());

    // configure the material data
    // http://ijbem.k.hosei.ac.jp/volume1/number1/pdf/ijbem_a4-10.pdf
    // http://www.ehow.com/about_6501091_conductivity-blood_.html
    //
    //  sigma_air = 0.0;
    //  sigma_bone = 0.018;
    //  sigma_brain = 0.25;
    //  sigma_plasma = 0.667;
    //
    // Layer 0 -> (sigma_bone, sigma_air)
    // Layer 1 -> (sigma_brain, sigma_bone)
    // Layer 2 -> (sigma_plasma, sigma_brain)
    //
    Assert(!parameters.material_data.empty(), ExcEmptyObject());
    material_data.read(parameters.material_data.c_str());

    // use Gauss-Legendre quadrature rule of order 4
    quadrature = QGauss<2>(4);

    // configure the solver control
    solver_control.log_frequency(1);
    solver_control.log_history(false);
    solver_control.log_result(true);
    solver_control.set_max_steps(100);
    solver_control.set_tolerance(1.0e-10);

    // read the boundary mesh for our domain
    Assert(!parameters.surface_mesh.empty(), ExcEmptyObject());
    sloc::read_ucd_mesh(parameters.surface_mesh.c_str(), tria);

    // check whether our other parameters are set
    Assert(!parameters.volume_mesh.empty(), ExcEmptyObject());
    Assert(!parameters.surface_phi.empty(), ExcEmptyObject());
    Assert(!parameters.volume_phi.empty(), ExcEmptyObject());

    // open logging stream
    if (!parameters.debug_logfile.empty())
        log.open(parameters.debug_logfile.c_str());

    // enumerate the basis functions, to figure out how many unknowns we've got
    dh.distribute_dofs(fe);

    // initialize vectors and matrices of our linear system
    const unsigned int n_dofs = dh.n_dofs();
    system_matrix.reinit(n_dofs, n_dofs);
    system_rhs.reinit(n_dofs);
    phi.reinit(n_dofs);
}

void BEM_ForwardProblem::assemble_system()
{
    using namespace std;

    deallog << "BEM_ForwardProblem::assemble_system() T=" << timer.wall_time() << std::endl;

    FEValues<2,3> fe_v(mapping, fe, quadrature,
                       update_values |
                       update_cell_normal_vectors |
                       update_quadrature_points |
                       update_JxW_values);

    const unsigned int n_dofs = dh.n_dofs();
    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
    dealii::Vector<double> local_matrix_row_i(fe.dofs_per_cell);

    typedef dealii::Point<3> Point3D;
    std::vector<Point3D> support_points(n_dofs);
    DoFTools::map_dofs_to_support_points<2,3>(mapping, dh, support_points);
    //sloc::write_points("support_points.dat", support_points);

    const bool debug = parameters.debug;

    // loop indices
    unsigned int i, j, q;

    // contribution from current dipole sources
    for (i = 0; i < n_dofs; ++i)
    {
        system_rhs(i) = 2 * dipole_sources.primary_contribution(support_points[i]);
    }

    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();


    if (debug)
    {
        log << "dh.n_dofs = " << dh.n_dofs() << endl;
        log << "fe.dofs_per_cell = " << fe.dofs_per_cell << endl;
        log << "fe_v.n_quadrature_points = " << fe_v.n_quadrature_points << endl;
    }

    ProgressTimer ptimer;
    cout << ptimer.header("cells");
    ptimer.start(tria.n_active_cells());
    for (unsigned int e = 0; cell != endc; ++cell, ++e)
    {
        fe_v.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const std::vector<Point3D>& q_points = fe_v.get_quadrature_points();
        const std::vector<Point3D>& normals = fe_v.get_cell_normal_vectors();

        const unsigned int mat_id = cell->material_id();
        const double sigma_int = material_data.get_sigma_int(mat_id);
        const double sigma_ext = material_data.get_sigma_ext(mat_id);
        const double sigma_avg = (sigma_int + sigma_ext) / 2;

        const double K =
            (1.0 / (4 * numbers::PI)) * ((sigma_int - sigma_ext) / sigma_avg);

        if (debug)
        {
            log << "------------\n";

            log << "K = " << K << endl;
            log << "cell "
                << cell->vertex_index(0) << " "
                << cell->vertex_index(1) << " "
                << cell->vertex_index(2) << " "
                << cell->vertex_index(3) << " "
                << endl;

            log << "local_dof_indices ";
            for (j = 0; j < fe.dofs_per_cell; ++j)
                log << local_dof_indices[j] << " ";
            log << endl;

            for (j = 0; j < fe.dofs_per_cell; ++j)
                log << "node " << cell->vertex_index(j) << " -> "
                    << cell->vertex(j)
                    << endl;

            if (debug && true)
            {
                log << "JxW q_points\n";
                for (q = 0; q < n_q_points; ++q)
                    log << fe_v.JxW(q) << ", " << q_points[q] << endl;
            }

            if (debug && true)
            {
                log << "shape_values\n";
                for (q = 0; q < n_q_points; ++q)
                {
                    for (j = 0; j < fe.dofs_per_cell; ++j)
                        log << fe_v.shape_value(j,q) << " ";
                    log << endl;
                }
            }

            if (debug && true)
            {
                log << "normals\n";
                for (q = 0; q < n_q_points; ++q)
                    log << normals[q] << endl;
            }

            if (debug && true)
            {
                log << "(R / R3)\n";
                for (q = 0; q < n_q_points; ++q)
                {
                    for (j = 0; j < fe.dofs_per_cell; ++j)
                    {
                        const Point<3> node_position = cell->vertex(j);
                        const Point<3> R = q_points[q] - node_position;
                        const double R3 = std::pow(R.square(), 1.5);
                        log << (R / R3) << ", ";
                    }
                    log << endl;
                }
            }
        }

        if (debug) log << "(R / R3) * normals\n";

        for (i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            DoFHandler<2,3>::active_cell_iterator cell2 = dh.begin_active();
            for (; cell2 != endc; ++cell2)
            {
                for (j = 0; j < fe.dofs_per_cell; ++j)
                {
                    if (i != cell2->vertex_index(j))
                        continue;

                    const Point<3> node_position = cell2->vertex(j);

                    if (debug) cout << "i=" << i << " j=" << j << "  ";

                    for (q = 0; q < n_q_points; ++q)
                    {
                        const Point<3> R = q_points[q] - node_position;

                        double R3 = std::pow(R.square(), 1.5);

                        double term = 
                            K * (R / R3) * normals[q] * fe_v.shape_value(j,q) * fe_v.JxW(q);

                        local_matrix_row_i(j) += term;

                        if (debug) cout << term << " + ";
                    }

                    if (debug) cout << "0 = " <<  local_matrix_row_i(j) << endl;
                }
            }

            for (j = 0; j < fe.dofs_per_cell; ++j)
            {
                system_matrix(i, local_dof_indices[j]) += -local_matrix_row_i(j);
            }
        }

        cout << ptimer.update(e);
    }
    cout << ptimer.update(tria.n_active_cells()) << endl;

    for (i = 0; i < n_dofs; ++i)
    {
        system_matrix(i,i) += 1.0;
    }

    if (debug)
    {
        sloc::write_matrix("system_matrix.dat", system_matrix);
        sloc::write_vector("system_rhs.dat", system_rhs);
    }
}

void BEM_ForwardProblem::solve_system()
{
    deallog << "BEM_ForwardProblem::solve_system() T=" << timer.wall_time() << std::endl;
    SolverGMRES<Vector<double> > solver(solver_control);
    solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity());
}

void BEM_ForwardProblem::output_results()
{
    deallog << "BEM_ForwardProblem::output_results() T=" << timer.wall_time() << std::endl;

    DataOut<2, DoFHandler<2,3> > data_out;
    data_out.attach_dof_handler(dh);
    data_out.add_data_vector(phi, "phi");
    data_out.build_patches(mapping, mapping.get_degree());

    std::stringstream ss;
    ss << parameters.surface_phi << ".vtk";
    std::string outfile = ss.str();
    std::ofstream os(outfile.c_str());
    data_out.write_vtk(os);
}

// ----------------------------------------------------------------------------

void BEM_ForwardProblem::compute_general_solution()
{
    deallog << "BEM_ForwardProblem::compute_general_solution() T="
            << timer.wall_time() << std::endl;

    Triangulation<3> g_tria;
    sloc::read_ucd_mesh(parameters.volume_mesh.c_str(), g_tria);

    FE_Q<3>         g_fe(1);
    DoFHandler<3>   g_dh(g_tria);
    Vector<double>  g_phi;

    g_dh.distribute_dofs(g_fe);

    g_phi.reinit(g_dh.n_dofs());

    FEValues<2,3> fe_v(mapping, fe, quadrature,
                       update_values |
                       update_cell_normal_vectors |
                       update_quadrature_points |
                       update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
    std::vector<double> local_phi(n_q_points);

    typedef dealii::Point<3> Point3D;
    std::vector<Point3D> g_support_points(g_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<3>(StaticMappingQ1<3>::mapping, g_dh, g_support_points);

    unsigned int i, q;

    // contribution from primary term
    for (i = 0; i < g_dh.n_dofs(); ++i)
    {
        g_phi(i) = 2 * dipole_sources.primary_contribution(g_support_points[i]);
    }

    //
    // contribution from surface integral term
    //
    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();

    for (; cell != endc; ++cell)
    {
        fe_v.reinit(cell);

        cell->get_dof_indices(local_dof_indices);
        fe_v.get_function_values(phi, local_phi);

        const std::vector<Point3D>& q_points = fe_v.get_quadrature_points();
        const std::vector<Point3D>& normals = fe_v.get_cell_normal_vectors();

        const unsigned int mat_id = cell->material_id();
        const double sigma_int = material_data.get_sigma_int(mat_id);
        const double sigma_ext = material_data.get_sigma_ext(mat_id);
        const double sigma_avg = (sigma_int + sigma_ext) / 2;

        for (q = 0; q < n_q_points; ++q)
        {
            const double K =
                (1.0 / (4 * numbers::PI)) * ((sigma_int - sigma_ext) / sigma_avg);

            // XXX: multiple dipole sources?
            const Point<3> dipole_source_position = dipole_sources(0).location;
            const Point<3> R = q_points[q] - dipole_source_position;
            const double R3 = std::pow(R.square(), 1.5);

            for (i = 0; i < g_dh.n_dofs(); ++i)
            {
                g_phi(i) += K * (local_phi[q] * (R / R3) * normals[q]) * fe_v.JxW(q);
            }
        }
    }

    DataOut<3> data_out;
    data_out.attach_dof_handler(g_dh);
    data_out.add_data_vector(g_phi, "g_phi");
    data_out.build_patches();

    std::stringstream ss;
    ss << parameters.volume_phi << ".vtk";
    std::string outfile = ss.str();
    std::ofstream os(outfile.c_str());
    data_out.write_vtk(os);
}

void BEM_ForwardProblem::compute_area()
{
    //
    // just a test.
    // we calculate the area by integrating the scalar 1 over the entire surface.
    //
    deallog << "BEM_ForwardProblem::compute_area() T=" << timer.wall_time() << std::endl;

    double area = 0;
    const double one = 1;

    FEValues<2,3> fe_v(mapping, fe, quadrature, update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);

    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();

    for (; cell != endc; ++cell)
    {
        fe_v.reinit(cell);
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int q = 0; q < n_q_points; ++q)
            area += one * fe_v.JxW(q);
    }

    deallog << "triangulation area = " << area << std::endl;
}

// ----------------------------------------------------------------------------
// EOF
