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
//#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
//#include <cstdlib>

/* other includes */
#include "io_dealii.h"

using namespace dealii;
using namespace sloc;

// ----------------------------------------------------------------------------

void BEM_ForwardProblem::Parameters::declare_parameters(ParameterHandler& prm)
{
    prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity");
    prm.declare_entry("dipole_sources", "", Patterns::Anything(), "Filename for current dipole sources data");
    prm.declare_entry("material_data", "", Patterns::Anything(), "Filename for material data");
    prm.declare_entry("surface_mesh", "", Patterns::Anything(), "Filename for surface mesh");
    prm.declare_entry("volume_mesh", "", Patterns::Anything(), "Filename for volume mesh");
}

void BEM_ForwardProblem::Parameters::get_parameters(ParameterHandler& prm)
{
    verbose = prm.get_bool("verbose");
    dipole_sources = prm.get("dipole_sources");
    material_data = prm.get("material_data");
    surface_mesh = prm.get("surface_mesh");
    volume_mesh = prm.get("volume_mesh");
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

    if (true)
    {
        cout << "dh.n_dofs = " << dh.n_dofs() << endl;
        cout << "fe.dofs_per_cell = " << fe.dofs_per_cell << endl;
        cout << "fe_v.n_quadrature_points = " << fe_v.n_quadrature_points << endl;
    }

    for (; cell != endc; ++cell)
    {
        // XXX: set up data dump for a single cell (for debugging)

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

        if (true)
        {
            cout << "------------\n";

            cout << "K = " << K << endl;
            cout << "cell " 
                 << cell->vertex_index(0) << " "
                 << cell->vertex_index(1) << " "
                 << cell->vertex_index(2) << " "
                 << cell->vertex_index(3) << " "
                 << endl;

            for (j = 0; j < fe.dofs_per_cell; ++j)
                cout << "node " << cell->vertex_index(j) << " -> "
                     << cell->vertex(j)
                     << endl;

            if (true)
            {
                cout << "JxW q_points\n";
                for (q = 0; q < n_q_points; ++q)
                    cout << fe_v.JxW(q) << ", " << q_points[q] << endl;
            }

            if (true)
            {
                cout << "shape_values\n";
                for (q = 0; q < n_q_points; ++q)
                {
                    for (j = 0; j < fe.dofs_per_cell; ++j)
                        cout << fe_v.shape_value(j,q) << " ";
                    cout << endl;
                }
            }

            if (true)
            {
                cout << "normals\n";
                for (q = 0; q < n_q_points; ++q)
                    cout << normals[q] << endl;
            }

            if (true)
            {
                cout << "(R / R3)\n";
                for (q = 0; q < n_q_points; ++q)
                {
                    for (j = 0; j < fe.dofs_per_cell; ++j)
                    {
                        const Point<3> node_position = cell->vertex(j);
                        const Point<3> R = q_points[q] - node_position;
                        const double R3 = std::pow(R.square(), 1.5);
                        cout << (R / R3) << ", ";
                    }
                    cout << endl;
                }
            }

        }

        for (i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            for (j = 0; j < fe.dofs_per_cell; ++j)
            {

                //if (i == cell->vertex_index(j)) continue;

                const Point<3> node_position = cell->vertex(j);

                for (q = 0; q < n_q_points; ++q)
                {
                    const Point<3> R = q_points[q] - node_position;

                    double R3 = std::pow(R.square(), 1.5);

                    double term = 
                        K * (R / R3) * normals[q] * fe_v.shape_value(j,q) * fe_v.JxW(q);

                    local_matrix_row_i(j) += term;

                    cout << std::setprecision(9) << term << " + ";
                }
                cout << "0 = " <<  local_matrix_row_i(j) << endl;
            }

            for (j = 0; j < fe.dofs_per_cell; ++j)
            {
                system_matrix(i, local_dof_indices[j]) += -local_matrix_row_i(j);
            }
        }
    }

    for (i = 0; i < n_dofs; ++i)
    {
        system_matrix(i,i) += 1.0;
    }

    if (true)
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

    std::ofstream outfile("phi.vtk");
    data_out.write_vtk(outfile);
}

// ----------------------------------------------------------------------------

void BEM_ForwardProblem::compute_general_solution()
{
    deallog << "BEM_ForwardProblem::compute_general_solution() T="
            << timer.wall_time() << std::endl;

    Triangulation<3> g_tria;
    sloc::read_ucd_mesh("tmp/head_volume.ucd", g_tria); // XXX: get filename from parameter file

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

    std::ofstream outfile("g_phi.vtk");
    data_out.write_vtk(outfile);
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
