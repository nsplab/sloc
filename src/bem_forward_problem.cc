#include "sloc/bem_forward_problem.h"

/* getfem includes */
#include <getfem/getfem_modeling.h>
#include <getfem/getfem_export.h>

/* deal.II includes */
// linear algebra
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
// various utilities
#include <deal.II/base/utilities.h>
#include <deal.II/base/table.h>

/* boost includes */
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

/* standard includes */
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <typeinfo>

/* our includes */
#include "sloc/io_dealii.h"
#include "sloc/progress_timer.h"

using namespace dealii;
using namespace sloc;
using namespace std;
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------
// XXX: move this elsewhere later
//

dealii::Point<3> cross(dealii::Point<3> A, dealii::Point<3> B)
{
    //
    //             | i    j   k |
    // C = A x B = | a0  a1  a2 |
    //             | b0  b1  b2 |
    //
    //   = (a1*b2 - a2*b1)i - (a0*b2 - a2*b0)j + (a0*b1 - a1*b0)k
    //
    double C[3];
    C[0] = +(A(1)*B(2) - A(2)*B(1));
    C[1] = -(A(0)*B(2) - A(2)*B(0));
    C[2] = +(A(0)*B(1) - A(1)*B(0));
    return dealii::Point<3>(C[0], C[1], C[2]);
} 

dealii::Point<3> triangle_normal(const bgeot::base_matrix& G)
{
    // columns of G are points in the corresponding simplex
    dealii::Point<3> A(G(0,0), G(1,0), G(2,0));
    dealii::Point<3> B(G(0,1), G(1,1), G(2,1));
    dealii::Point<3> C(G(0,2), G(1,2), G(2,2));
    dealii::Point<3> N = cross(C-A, B-A);
    return N / std::sqrt(N.square());
}

template <class CONT>
void copy_to_cout(const CONT& c, const char *delim)
{
    std::copy(c.begin(), c.end(),
        std::ostream_iterator<typename CONT::value_type>(std::cout, delim));
}

// ----------------------------------------------------------------------------

BEM_Forward_Problem::BEM_Forward_Problem(const Parameters& parameters)
    : parameters(parameters)    // parameters class
    , mf(surface_mesh)          // finite element method on mesh
    , mim(surface_mesh)         // integration method on mesh
    , verbose(true)             // print progress by default
    , debug(false)              // don't print out debug messages by default
{
    timer.start();
}

BEM_Forward_Problem::~BEM_Forward_Problem()
{
}

void BEM_Forward_Problem::run()
{
    deallog << "BEM_Forward_Problem::run() T=" << timer.wall_time() << endl;
    configure();
    compute_area();
    assemble_system();
    solve_system();
    output_results();
    //compute_general_solution();
    deallog << "last timer = " << timer.wall_time() << endl;
}

// ----------------------------------------------------------------------------

void BEM_Forward_Problem::Parameters::declare_parameters(ParameterHandler& prm)
{
    prm.declare_entry("surface_mesh", "", Patterns::Anything(), "Filename for surface mesh");
    prm.declare_entry("surface_mesh_materials", "", Patterns::Anything(), "Filename for surface mesh materials");
    prm.declare_entry("material_data", "", Patterns::Anything(), "Filename for material data");
    prm.declare_entry("dipole_sources", "", Patterns::Anything(), "Filename for current dipole sources data");

    prm.declare_entry("output_vtk", "", Patterns::Anything(), "XXX");
    prm.declare_entry("output_phi", "", Patterns::Anything(), "XXX");

    prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level");
    prm.declare_entry("debug", "false", Patterns::Bool(), "Output debug information");
    prm.declare_entry("logfile", "debug.log", Patterns::Anything(), "Filename for debugging logfile");
}

void BEM_Forward_Problem::Parameters::get_parameters(ParameterHandler& prm)
{
    surface_mesh = prm.get("surface_mesh");
    surface_mesh_materials = prm.get("surface_mesh_materials");
    material_data = prm.get("material_data");
    dipole_sources = prm.get("dipole_sources");

    output_vtk = prm.get("output_vtk");
    output_phi = prm.get("output_phi");

    verbose = prm.get_bool("verbose");
    debug = prm.get_bool("debug");
    logfile = prm.get("logfile");
}

// ----------------------------------------------------------------------------

void BEM_Forward_Problem::configure()
{
    deallog << "BEM_Forward_Problem::configure() T=" << timer.wall_time() << endl;

    // keep a stringstream handy so we can build error messages before throwing
    stringstream ss;

    // remember these two flags
    debug = parameters.debug;
    verbose = parameters.verbose;

    // configure the dipole sources (but only if they haven't already been configured)
    if (dipole_sources.n_sources() == 0)
    {
        Assert(!parameters.dipole_sources.empty(), ExcEmptyObject());
        if (fs::exists(parameters.dipole_sources))
        {
            if (debug) cerr << "Reading " << parameters.dipole_sources << endl;
            dipole_sources.read(parameters.dipole_sources.c_str());
        }
        else
        {
            ss << "Invalid file " << parameters.dipole_sources;
            throw std::invalid_argument(ss.str());
        }
    }

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
    if (fs::exists(parameters.material_data))
    {
        if (debug) cerr << "Reading " << parameters.material_data << endl;
        material_data.read_sigma(parameters.material_data.c_str());
    }
    else
    {
        ss << "Invalid file " << parameters.material_data;
        throw std::invalid_argument(ss.str());
    }

    // initialize finite element space
    pf = getfem::fem_descriptor("FEM_PK(2,1)");

    // initialize geometric transformation
    pgt = bgeot::simplex_geotrans(2,1);

    // initialize integration rule
    pim = getfem::classical_approx_im(pgt, 6);

    // read the surface mesh we'll be integrating over
    Assert(!parameters.surface_mesh.empty(), ExcEmptyObject());
    if (fs::exists(parameters.surface_mesh))
    {
        if (debug) cerr << "Reading " << parameters.surface_mesh << endl;
        surface_mesh.read_from_file(parameters.surface_mesh);

        if (debug && verbose)
        {
            cerr << "  found " << surface_mesh.nb_points() << " points, and "
                 << surface_mesh.nb_convex() << " cells" << endl;
        }
    }
    else
    {
        ss << "Invalid file " << parameters.surface_mesh;
        throw std::invalid_argument(ss.str());
    }

    // assign the finite element
    mf.set_finite_element(pf);
    mim.set_integration_method(pim);

    // read the surface mesh materials
    Assert(!parameters.surface_mesh_materials.empty(), ExcEmptyObject());
    if (fs::exists(parameters.surface_mesh_materials))
    {
        if (debug) cerr << "Reading " << parameters.surface_mesh_materials << endl;
        material_data.read_materials(parameters.surface_mesh_materials.c_str());
    }
    else
    {
        ss << "Invalid file " << parameters.surface_mesh_materials;
        throw std::invalid_argument(ss.str());
    }

    // configure the solver control
    solver_control.log_frequency(1);
    solver_control.log_history(false);
    solver_control.log_result(true);
    solver_control.set_max_steps(100);
    solver_control.set_tolerance(1.0e-10);

    // open logging stream, but only in debug mode
    if (debug && !parameters.logfile.empty())
        log.open(parameters.logfile.c_str());

    // initialize vectors and matrices of our linear system
    const unsigned int n_dofs = mf.nb_basic_dof();
    system_matrix.reinit(n_dofs, n_dofs);
    system_rhs.reinit(n_dofs);
    phi.reinit(n_dofs);
}


void BEM_Forward_Problem::assemble_system()
{
    deallog << "BEM_Forward_Problem::assemble_system() T=" << timer.wall_time() << endl;
    
    const bool debug_assembly = true;

    const unsigned int n_dofs = mf.nb_basic_dof();
    const int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
    assert(fe_dofs_per_cell == 3);

    std::vector<unsigned int> local_dof_indices(fe_dofs_per_cell);
    dealii::Vector<double> local_matrix_row_i(fe_dofs_per_cell);

    // loop indices
    unsigned int i, j, q;

    // for holding the local points
    bgeot::base_matrix G;

    // zero out the system matrix, so we can call this method repeatedly
    system_matrix.reset_values();

    if (debug)
    {
        log << "n_dofs = " << n_dofs << endl;
        log << "fe_dofs_per_cell = " << fe_dofs_per_cell << endl;
    }


    //
    // contribution from current dipole sources
    //
    for (i = 0; i < n_dofs; ++i)
    {
        bgeot::base_node p = mf.point_of_basic_dof(i);
        dealii::Point<3> pt(p[0], p[1], p[2]);
        system_rhs(i) = 2 * dipole_sources.primary_contribution(pt);
    }


    ProgressTimer ptimer;
    long e = 0;
    if (verbose)
    {
        cerr << ptimer.header("cells");
        ptimer.start(surface_mesh.nb_convex());
        //log << "surface_mesh.nb_convex() = " << surface_mesh.nb_convex() << endl;
    }

    //
    // contribution from surface integral terms
    //
    for (dal::bv_visitor cv(surface_mesh.convex_index()); !cv.finished(); ++cv)
    {

        bgeot::pgeometric_trans pgt = surface_mesh.trans_of_convex(cv);
        getfem::papprox_integration pai = getfem::get_approx_im_or_fail(mim.int_method_of_element(cv));
        getfem::pfem pf = mf.fem_of_element(cv);

        bgeot::vectors_to_base_matrix(G, surface_mesh.points_of_convex(cv));

        // for curved elements, this might depend on q. for our linear triangle,
        // it's constant so we factor this term outside of the inner loops below
        dealii::Point<3> normal_q = triangle_normal(G);

        getfem::fem_interpolation_context ctx(pgt, pf, bgeot::base_node(3), G, cv);

        // get current element's dof indices
        getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);

        // retrieve material data
        const unsigned int mat_id = material_data.get_material_id(cv);
        const double sigma_int = material_data.get_sigma_int(mat_id);
        const double sigma_ext = material_data.get_sigma_ext(mat_id);
        const double sigma_avg = (sigma_int + sigma_ext) / 2;

        const double K = (1.0 / (4 * numbers::PI)) * (sigma_int - sigma_ext) / sigma_avg;

        if (debug_assembly)
        {
            log << "convex " << cv << endl;
            //log << "convex " << cv << " - " << typeid(cv).name() << endl;

            log << "G = " << G;
            log << "normal_q = " << normal_q << endl;
            log << "mat_id = " << mat_id << endl;
            log << "sigma_int = " << sigma_int << endl;
            log << "sigma_ext = " << sigma_ext << endl;
            log << "sigma_avg = " << sigma_avg << endl;

            //copy_to_cout(elt_dof_indices, " "); cout << endl;
            log << "elt_dof_indices = ";
            std::copy(elt_dof_indices.begin(), elt_dof_indices.end(),
                std::ostream_iterator<getfem::mesh_fem::ind_dof_ct::value_type>(log, " "));
            log << endl;
        }

        for (i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            bgeot::base_node p = mf.point_of_basic_dof(i);
            const Point<3> node_position(p[0], p[1], p[2]);
            if (debug_assembly)
            {
                log << "  dof " << i << ", pt [" << node_position << "], nq = " << pai->nb_points_on_convex() << endl;
            }

            for (q = 0; q < pai->nb_points_on_convex(); ++q)
            {
                // set current point (called xref) on the reference element
                ctx.set_xref(pai->point(q));

                // evaluate shape functions at xref
                getfem::base_tensor t;
                ctx.base_value(t);

                // calculate x (geometric transform applied to xref)
                bgeot::base_node x;
                x = pgt->transform(pai->point(q), G);

                const Point<3> q_point(x[0], x[1], x[2]);
                const Point<3> R = q_point - node_position;
                const double R3 = std::pow(R.square(), 1.5);
                const double JxW = pai->coeff(q) * ctx.J();

                for (j = 0; j < fe_dofs_per_cell; ++j)
                {
                    // assemble the matrix!
                    const double shape_value_j = t(j,0);
                    const double term = K * (R / R3) * normal_q * shape_value_j * JxW;
                    local_matrix_row_i(j) += term;
                }

                if (debug_assembly)
                {
                    log << "    q = " << q << endl;
                    log << "    xref = " << pai->point(q) << endl;
                    log << "    N(xref) = [" << t(0,0) << ", " << t(1,0) << ", " << t(2,0) << "]" << endl;
                    log << "    x = " << x << endl;
                    log << "    J = " << ctx.J() << endl;
                    log << "    R = " << R << endl;
                    log << "    R3 = " << R3 << endl;
                    //log << "    local_matrix_row_" << i << " = " << local_matrix_row_i << endl;
                    for (j = 0; j < fe_dofs_per_cell; ++j)
                        log << "    local_matrix(" << i << "," << j << ") = " << local_matrix_row_i(j) << endl;
                    log << "    ------------------------" << endl;
                }
            }
            if (debug_assembly)
            {
                for (j = 0; j < fe_dofs_per_cell; ++j)
                    log << "    local_matrix(" << i << "," << j << ") = " << local_matrix_row_i(j) << endl;
            }

            for (j = 0; j < fe_dofs_per_cell; ++j)
            {
                system_matrix(i, elt_dof_indices[j]) += -local_matrix_row_i(j);
            }
        }

        if (verbose) cerr << ptimer.update(e++);
    }
    if (verbose) cerr << ptimer.update(e) << endl;

    // Finally, add the identity matrix to system_matrix
    for (i = 0; i < n_dofs; ++i)
        system_matrix(i,i) += 1.0;

    // Write the linear system to disk, for debugging
    if (debug)
    {
        sloc::write_vector("system_rhs.dat", system_rhs);

        if (n_dofs < 1000)
            sloc::write_matrix("system_matrix.dat", system_matrix);
    }

}

void BEM_Forward_Problem::solve_system()
{
    //
    // See http://en.wikipedia.org/wiki/Generalized_minimal_residual_method
    //
    deallog << "BEM_Forward_Problem::solve_system() T=" << timer.wall_time() << endl;
    SolverGMRES<Vector<double> > solver(solver_control);
    solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity());
}

void BEM_Forward_Problem::output_results()
{
    deallog << "BEM_Forward_Problem::output_results() T=" << timer.wall_time() << endl;

    if (!parameters.output_phi.empty())
    {
        sloc::write_vector(parameters.output_phi.c_str(), phi);
    }

    // XXX: what is this for, again?

    if (!parameters.output_vtk.empty())
    {
        // copy the solution vector phi into a getfem vector that we can save
        getfem::modeling_standard_plain_vector phi_data;
        std::copy(phi.begin(), phi.end(), std::back_inserter(phi_data));
        //mf.write_to_file("tetra.meshfem", true);

        // write out vtk file
        getfem::vtk_export exp(parameters.output_vtk.c_str(), true);
        exp.exporting(mf);
        exp.write_point_data(mf, phi_data, "phi");
    }

}

void BEM_Forward_Problem::compute_general_solution()
{
    deallog << "BEM_Forward_Problem::compute_general_solution() T=" << timer.wall_time() << endl;
    Assert(false, ExcNotImplemented());
}

void BEM_Forward_Problem::compute_area()
{
    //
    // just a test.
    // we calculate the area by integrating the scalar 1 over the entire surface
    //
    deallog << "BEM_Forward_Problem::compute_area() T=" << timer.wall_time() << endl;

    double area = 0;
    const double one = 1.0;

    bgeot::base_matrix G;

    for (dal::bv_visitor cv(surface_mesh.convex_index()); !cv.finished(); ++cv)
    {
        bgeot::pgeometric_trans pgt = surface_mesh.trans_of_convex(cv);
        getfem::papprox_integration pai = getfem::get_approx_im_or_fail(mim.int_method_of_element(cv));
        getfem::pfem pf = mf.fem_of_element(cv);

        bgeot::vectors_to_base_matrix(G, surface_mesh.points_of_convex(cv));

        getfem::fem_interpolation_context ctx(pgt, pf, bgeot::base_node(3), G, cv);

        for (unsigned int q = 0; q < pai->nb_points_on_convex(); ++q)
        {
            ctx.set_xref(pai->point(q));
            area += one * pai->coeff(q) * ctx.J();
        }
    }

    if (debug) cerr << "surface mesh area = " << area << endl;
}

// ----------------------------------------------------------------------------
// EOF
