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
#include "sloc/geometry.h"

using namespace dealii;
using namespace sloc;
using namespace std;
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------
// XXX: move this elsewhere later
//

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
    //solver_control.set_tolerance(1.0e-10);
    solver_control.set_tolerance(6.0e-5);

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
    const unsigned int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
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
    ofstream mapping("mapping.dat");
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
        getfem::mesh_fem::ind_dof_face_ct elt_dof_f_indices = mf.ind_basic_dof_of_face_of_element(cv,0);
        //getfem::mesh cv_mesh = mf.linked_mesh();
	//getfem::mesh::ref_mesh_pt_ct mesh_ptc_idx = cv_mesh.points_of_convex(cv);
	std::vector< getfem::size_type > global_ind;
	mf.get_global_dof_index(global_ind);
	for (size_t tidx=0; tidx<global_ind.size(); tidx++)
		cout<<"tidx: "<<tidx<<" "<<global_ind[tidx]<<endl;

        // retrieve material data
        const unsigned int mat_id = material_data.get_material_id(cv);
        const double sigma_int = material_data.get_sigma_int(mat_id);
        const double sigma_ext = material_data.get_sigma_ext(mat_id);
        const double sigma_avg = (sigma_int + sigma_ext) / 2.0;

        const double K = (1.0 / (4 * numbers::PI)) * (sigma_int - sigma_ext) / sigma_avg;


        //if (debug_assembly)
        {
            cout << "convex " << cv << endl;
            //log << "convex " << cv << " - " << typeid(cv).name() << endl;

            cout << "G = " << G;
            cout << "normal_q = " << normal_q << endl;
            cout << "mat_id = " << mat_id << endl;
            cout << "sigma_int = " << sigma_int << endl;
            cout << "sigma_ext = " << sigma_ext << endl;
            cout << "sigma_avg = " << sigma_avg << endl;

            //copy_to_cout(elt_dof_indices, " "); cout << endl;
            log << "elt_dof_indices = ";
            std::copy(elt_dof_indices.begin(), elt_dof_indices.end(),
                std::ostream_iterator<getfem::mesh_fem::ind_dof_ct::value_type>(log, " "));
            log << endl;
        }

    //
    // contribution from current dipole sources
    //
	mapping<<cv;
    for (j = 0; j < fe_dofs_per_cell; ++j)
    {
        bgeot::base_node p = mf.point_of_basic_dof(elt_dof_indices[j]);
        dealii::Point<3> pt(p[0], p[1], p[2]);
	cout<<"cv: "<<cv<<endl;
	cout<<"j: "<<j<<endl;
	cout<<"elt dof indx: "<<elt_dof_indices[j]<<endl;
	cout<<"elt dof face indx: "<<elt_dof_f_indices[j]<<endl;
	cout<<"elt dof mesh: "<<global_ind[elt_dof_indices[j]]<<endl;
	cout<<"loc: "<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
	cout<<"red: "<<mf.is_reduced()<<endl;
        system_rhs(elt_dof_indices[j]) = dipole_sources.primary_contribution(pt) / sigma_avg;
	mapping<<" "<<elt_dof_indices[j];
    }
	mapping<<endl;

        if (verbose) cout<<"n_dofs: "<<n_dofs<<endl;
        if (verbose) cout<<"pai->nb_points_on_convex(): "<<(pai->nb_points_on_convex())<<endl;
        if (verbose) cout<<"fe_dofs_per_cell: "<<fe_dofs_per_cell<<endl;
        for (i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            //
            // contribution from current dipole sources
            //
            bgeot::base_node p = mf.point_of_basic_dof(i);
            //dealii::Point<3> pt(p[0], p[1], p[2]);
            //cout<<"sigma_avg: "<<sigma_avg<<endl;
            //system_rhs(i) = dipole_sources.primary_contribution(pt) / sigma_avg;
            //cout<<"rhs "<<i<<": "<<dipole_sources.primary_contribution(pt)<<endl;
            //cout<<"rhs "<<i<<": "<<system_rhs(i)<<endl;


            const Point<3> node_position(p[0], p[1], p[2]);
	cout<<"i: "<<i<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
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
                system_matrix(i, elt_dof_indices[j]) += local_matrix_row_i(j);
                //cout << "  k=0 cv=" << cv << " i=" << i << " j=" << j << " M(" << i << "," << elt_dof_indices[j] << ") += " << -local_matrix_row_i(j) << endl;
            }
        }

        if (verbose) cerr << ptimer.update(e++);
    }
    if (verbose) cerr << ptimer.update(e) << endl;

    // Finally, add the identity matrix to system_matrix
    for (i = 0; i < n_dofs; ++i)
        system_matrix(i,i) += 1.0;

    // Write the linear system to disk, for debugging
    //if (debug)
    {
        sloc::write_vector("system_rhs_s.dat", system_rhs);

        if (n_dofs < 1000)
            sloc::write_matrix("system_matrix_s.dat", system_matrix);
    }

}

void BEM_Forward_Problem::solve_system()
{
    //
    // See http://en.wikipedia.org/wiki/Generalized_minimal_residual_method
    //
    deallog << "BEM_Forward_Problem::solve_system() T=" << timer.wall_time() << endl;
    try
    {
        cout<<"before SolverGMRES"<<endl;
        SolverGMRES<Vector<double> > solver(solver_control);
        cout<<"after SolverGMRES"<<endl;
        solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity());
        cout<<"after solve"<<endl;
    }
    catch (std::exception& exc)
    {
        //cout<<exc.what()<<endl;
        ofstream out;
        out.open("head.threshold_error", std::ofstream::out | std::ofstream::app);
        out<<exc.what()<<endl;
        out.close();
        //sloc::write_matrix("system_matrix.dat", system_matrix);
        //sloc::write_vector("system_rhs.dat", system_rhs);
        //throw exc;
    }
}

void BEM_Forward_Problem::output_results()
{
    deallog << "BEM_Forward_Problem::output_results() T=" << timer.wall_time() << endl;

    if (!parameters.output_phi.empty())
    {
        sloc::write_vector(parameters.output_phi.c_str(), phi);
    }

    if (!parameters.output_vtk.empty())
    {
        // copy the solution vector phi into a getfem vector that we can save
        getfem::modeling_standard_plain_vector phi_data;
        std::copy(phi.begin(), phi.end(), std::back_inserter(phi_data));

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
