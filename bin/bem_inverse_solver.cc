// bem_inverse_solver.cc
// ----------------------------------------------------------------------------

/* Inverse Solver
 *
 *  Inputs are:
 *  (1) parameter file for running forward problem
 *      - surface mesh
 *      - verbose
 *  (2) an electrode file containing the node ids and potential values
 *      at the electrode locations
 *
 *  Outputs are:
 *  - a sources.dipole file containing a list of inferred dipoles
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include <sloc/dipole_sources.h>
#include <sloc/simplex_search.h>
#include <sloc/bem_forward_problem.h>
//#include <sloc/bem_inverse_problem.h>
#include <sloc/io_dealii.h>

namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------

class ForwardProblem : public sloc::BEM_Forward_Problem
{

public:

    ForwardProblem(const sloc::BEM_Forward_Problem::Parameters parameters)
        : sloc::BEM_Forward_Problem(parameters)
    {
    }

    void configure()
    {
        // first, create a dummy dipole source
        dealii::Point<3> dummy_x(0,0,0), dummy_p(0,0,1);
        dipole_sources.add_source(dummy_x, dummy_p);

        // configure the rest of the parameters
        sloc::BEM_Forward_Problem::configure();
    }

    void solve(double x, double y, double z, double px, double py, double pz)
    {
        // Assume we only have a single dipole
        Assert(dipole_sources.n_sources() == 1, dealii::ExcInvalidState());
        sloc::DipoleSource *src = dipole_sources._sources[0];
        src->location = dealii::Point<3>(x,y,z);
        src->dipole = dealii::Point<3>(px,py,pz);

        // build the linear system, and solve it
        this->assemble_system();
        this->solve_system();
    }

    void solve(double *points, double *dipoles)
    {
        for (unsigned int i = 0; i < dipole_sources.n_sources(); i++)
        {
            sloc::DipoleSource *src = dipole_sources._sources[i];
            src->location = dealii::Point<3>(points[3*i+0], points[3*i+1], points[3*i+2]);
            src->dipole = dealii::Point<3>(dipoles[3*i+0], dipoles[3*i+1], dipoles[3*i+2]);
        }
        this->assemble_system();
        this->solve_system();
    }

    double get_phi(int k)
    {
        return phi(k);
    }

};


// ----------------------------------------------------------------------------

class InverseProblem
{
public:
    class Parameters
    {
    public:
        static void declare_parameters(dealii::ParameterHandler& prm)
        {
            using namespace dealii;

            prm.declare_entry("electrodes", "electrodes.dat", Patterns::Anything(), "Input file with electrode potential measurements");
            prm.declare_entry("output_sources", "sources.dipole", Patterns::Anything(), "Final output of source localization");
            prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level of inverse solver");
            prm.declare_entry("debug", "false", Patterns::Bool(), "Debugging information for inverse solver");

            prm.enter_subsection("Forward Problem Parameters");
            {
                prm.declare_entry("surface_mesh", "", Patterns::Anything(), "Filename for surface mesh");
                prm.declare_entry("surface_mesh_materials", "", Patterns::Anything(), "Filename for surface mesh materials");
                prm.declare_entry("material_data", "", Patterns::Anything(), "Filename for material data");
                prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level of forward solver");
                prm.declare_entry("debug", "false", Patterns::Bool(), "Debugging information for forward solver");
            }
            prm.leave_subsection();

            prm.enter_subsection("Simplex Search Parameters");
            {
                prm.declare_entry("initial_search_point", "0,0,0", Patterns::List(Patterns::Double(), 3, 3), "Initial search point for simplex search");
                prm.declare_entry("initial_search_radius", "0.1", Patterns::Double(), "Initial search radius for simplex search");
                prm.declare_entry("tolerance", "1e-8", Patterns::Double(), "Convergence criterion for simplex search");
                prm.declare_entry("max_iterations", "1000", Patterns::Integer(), "Upper limit on number of iterations for simplex search");
                prm.declare_entry("reflection_coefficient", "1.0", Patterns::Double(), "");
                prm.declare_entry("contraction_coefficient", "0.5", Patterns::Double(), "");
                prm.declare_entry("expansion_coefficient", "2.0", Patterns::Double(), "");
                prm.declare_entry("reduction_coefficient", "0.5", Patterns::Double(), "");
                prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level for simplex search");
                prm.declare_entry("debug", "false", Patterns::Bool(), "Debugging information for simplex search");
            }
            prm.leave_subsection();
        }

        void get_parameters(dealii::ParameterHandler& prm)
        {
            electrodes = prm.get("electrodes");
            output_sources = prm.get("output_sources");
            verbose = prm.get_bool("verbose");
            debug = prm.get_bool("debug");

            prm.enter_subsection("Forward Problem Parameters");
            {
                surface_mesh = prm.get("surface_mesh");
                surface_mesh_materials = prm.get("surface_mesh_materials");
                material_data = prm.get("material_data");
                fwd_verbose = prm.get_bool("verbose");
                fwd_debug = prm.get_bool("debug");
            }
            prm.leave_subsection();

            prm.enter_subsection("Simplex Search Parameters");
            {
                // convert string with comma separated values into a vector<double>
                using namespace std;
                using namespace boost::algorithm;
                vector<string> coords;
                std::string str = prm.get("initial_search_point");
                split(coords, str, is_any_of(", "), token_compress_on);
                for (vector<string>::iterator it = coords.begin(); it != coords.end(); ++it)
                    initial_search_point.push_back(boost::lexical_cast<double>(*it));

                // read the rest of the parameters
                initial_search_radius = prm.get_double("initial_search_radius");
                tolerance = prm.get_double("tolerance");
                max_iterations = prm.get_integer("max_iterations");
                alpha = prm.get_double("reflection_coefficient");
                beta = prm.get_double("contraction_coefficient");
                gamma = prm.get_double("expansion_coefficient");
                sigma = prm.get_double("reduction_coefficient");
                simplex_verbose = prm.get_bool("verbose");
                simplex_debug = prm.get_bool("debug");
            }
            prm.leave_subsection();
        }

        void print() const
        {
            using namespace std;
            cout << "electrodes = " << electrodes << endl;
            cout << "output_sources = " << output_sources << endl;
            cout << "surface_mesh = " << surface_mesh << endl;
            cout << "surface_mesh_materials = " << surface_mesh_materials << endl;
            cout << "material_data = " << material_data << endl;
            cout << "initial_search_point = " << initial_search_point[0] << ", " << initial_search_point[1] << ", " << initial_search_point[2] << endl;
            cout << "initial_search_radius = " << initial_search_radius << endl;
            cout << "alpha = " << alpha << endl;
            cout << "beta = " << beta << endl;
            cout << "gamma = " << gamma << endl;
            cout << "sigma = " << sigma << endl;
            cout << "max_iterations = " << max_iterations << endl;
            cout << "verbose = " << verbose << endl;
            cout << "fwd_verbose = " << fwd_verbose << endl;
            cout << "simplex_verbose = " << simplex_verbose << endl;
            cout << "debug = " << debug << endl;
            cout << "fwd_debug = " << fwd_debug << endl;
            cout << "simplex_debug = " << simplex_debug << endl;
        }

    public:
        std::string electrodes;
        std::string output_sources;
        std::string surface_mesh;
        std::string surface_mesh_materials;
        std::string material_data;
        std::vector<double> initial_search_point;
        double initial_search_radius;
        double alpha, beta, gamma, sigma;
        unsigned int max_iterations;
        double tolerance;
        bool verbose, fwd_verbose, simplex_verbose;
        bool debug, fwd_debug, simplex_debug;
    };

public:

    InverseProblem(Parameters &parameters) : parameters(parameters)
    {
        // disable dealii logging (fwd solver needs to be quiet)
        dealii::deallog.depth_console(0);

        // a few sanity checks
        Assert(!parameters.electrodes.empty(), dealii::ExcEmptyObject());
        Assert(!parameters.output_sources.empty(), dealii::ExcEmptyObject());
        Assert(!parameters.surface_mesh.empty(), dealii::ExcEmptyObject());
        Assert(!parameters.surface_mesh_materials.empty(), dealii::ExcEmptyObject());
        Assert(!parameters.material_data.empty(), dealii::ExcEmptyObject());

        // Check that electrodes file exists (XXX: throwing exception in constructor?? fix this)
        if (!fs::exists(parameters.electrodes))
        {
            std::stringstream ss;
            ss << "Invalid file " << parameters.electrodes << std::endl;
            throw std::invalid_argument(ss.str());
        }

        //
        // Initialize electrodes (XXX: refactor this section of code)
        //

        std::ifstream is;
        is.open(parameters.electrodes.c_str()); // XXX: check that open succeeded

        is >> num_electrodes;
        const int M = num_electrodes;

        phi_measured.reinit(M);
        electrode_dof_index.reserve(M);
        for (int m = 0; m < M; m++)
        {
            unsigned int dof_index;
            double phi_value;

            is >> dof_index;
            is >> phi_value;

            phi_measured(m) = phi_value;
            electrode_dof_index.push_back(dof_index);
        }

        //
        // Initialize lead-field matrix and related vectors
        //

        L.reinit(M, 3);
        L_pinv.reinit(3,M);
        dipole.reinit(3);
    }

    double cost_for_single_source_at(double pt[3], double dp[3])
    {
        sloc::BEM_Forward_Problem::Parameters fwd_prm;
        fwd_prm.surface_mesh = parameters.surface_mesh;
        fwd_prm.surface_mesh_materials = parameters.surface_mesh_materials;
        fwd_prm.material_data = parameters.material_data;
        fwd_prm.dipole_sources = "";
        fwd_prm.verbose = parameters.fwd_verbose;
        fwd_prm.debug = parameters.fwd_debug;
        fwd_prm.logfile = "";

        ForwardProblem forward_problem(fwd_prm);
        forward_problem.configure();

        int m;
        const int M = num_electrodes;

        double x = pt[0];
        double y = pt[1];
        double z = pt[2];

        // find L[:,0]
        forward_problem.solve(x, y, z, 1, 0, 0);
        for (m = 0; m < M; m++)
            L(m,0) = forward_problem.get_phi(electrode_dof_index[m]);

        // find L[:,1]
        forward_problem.solve(x, y, z, 0, 1, 0);
        for (m = 0; m < M; m++)
            L(m,1) = forward_problem.get_phi(electrode_dof_index[m]);

        // find L[:,2]
        forward_problem.solve(x, y, z, 0, 0, 1);
        for (m = 0; m < M; m++)
            L(m,2) = forward_problem.get_phi(electrode_dof_index[m]);

        // find the pseudo-inverse of L, and store it in the matrix L_pinv
        L_pinv.left_invert(L);

        // calculate: dipole = pinv(L) * phi_measured
        L_pinv.vmult(dipole, phi_measured);
        dp[0] = dipole(0);
        dp[1] = dipole(1);
        dp[2] = dipole(2);

        // compute cost ||phi - L * dipole||
        double cost = 0.0;
        for (m = 0; m < M; m++)
        {
            double dphi = phi_measured(m) - dipole(0) * L(m,0)
                                          - dipole(1) * L(m,1)
                                          - dipole(2) * L(m,2);
            cost += dphi * dphi;
        }
        cost = sqrt(cost);

        return cost;
    }

    double cost_for_sources_at(double * /*positions*/)
    {
        // XXX: implement this method
        double cost = 0.0;
        return cost ;
    }


public:

    const Parameters& parameters;

    unsigned int num_electrodes;
    std::vector<unsigned int> electrode_dof_index;
    dealii::Vector<double> phi_measured;

    dealii::FullMatrix<double> L;
    dealii::FullMatrix<double> L_pinv;
    dealii::Vector<double> dipole;

};

/* global objects! */
static InverseProblem *g_inverse_problem = 0;
static InverseProblem::Parameters *g_parameters = 0;

// ----------------------------------------------------------------------------

static void init_globals(char *parameters_file)
{
    g_parameters = new InverseProblem::Parameters;
    dealii::ParameterHandler parameter_handler;
    g_parameters->declare_parameters(parameter_handler);
    if (!parameter_handler.read_input(parameters_file)) exit(1);
    g_parameters->get_parameters(parameter_handler);
    g_inverse_problem = new InverseProblem(*g_parameters);
}

static void cleanup_globals()
{
    if (g_parameters)
    {
        delete g_parameters;
        g_parameters = 0;
    }

    if (g_inverse_problem)
    {
        delete g_inverse_problem;
        g_inverse_problem = 0;
    }
}

double cost_function(double x[3])
{
    double p[3];
    return g_inverse_problem->cost_for_single_source_at(x,p);
}

void minimize_cost_function()
{
    const InverseProblem::Parameters& prm = g_inverse_problem->parameters;

    const int dim = 3;
    sloc::NelderMead::SimplexSearch<dim> simplex_search;
    simplex_search.F = cost_function;
    simplex_search.tol = prm.tolerance;
    simplex_search.max_iterations = prm.max_iterations;
    simplex_search.alpha = prm.alpha;
    simplex_search.beta = prm.beta;
    simplex_search.gamma = prm.gamma;
    simplex_search.sigma = prm.sigma;
    simplex_search.debug = prm.simplex_debug;
    simplex_search.verbose = prm.simplex_verbose;

    double neighborhood = prm.initial_search_radius;
    double initial_point[dim] = {
        prm.initial_search_point[0],
        prm.initial_search_point[1],
        prm.initial_search_point[2]
    };
    double initial_simplex[4][dim];
    sloc::NelderMead::SimplexSearch<dim>::make_wedge(neighborhood, initial_point, initial_simplex);

    double final_cost;
    double final_point[dim];
    int ret = simplex_search.run(initial_simplex, final_point, final_cost);

    // XXX: can we get final_dipole without evaluating cost function again?
    double final_dipole[dim];
    g_inverse_problem->cost_for_single_source_at(final_point, final_dipole);

    if (ret == 0)
    {
        std::cout << "Found minimum at ("
                  << final_point[0] << ", "
                  << final_point[1] << ", "
                  << final_point[2] << ") with cost "
                  << final_cost
                  << std::endl;
    }

    sloc::DipoleSources sources;
    dealii::Point<3> loc(final_point[0], final_point[1], final_point[2]);
    dealii::Point<3> dip(final_dipole[0], final_dipole[1], final_dipole[2]);
    sources.add_source(loc, dip);
    sources.write(prm.output_sources.c_str());

    // done!
    return;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    using namespace std;

    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <input.prm>" << endl;
        return 1;
    }

    if (!fs::exists(argv[1]))
    {
        cerr << "Error: parameter file '" << argv[1] << "' does not exist!" << endl;
        return 2;
    }

    string sep = "----------------------------------------------------";

    try
    {
        cerr << "bem_inverse_solver" << endl;

        // XXX: use this one later
        //BEM_Inverse_Problem inverse_problem;
        //inverse_problem.run();

        // XXX: for now, do things the old way
        init_globals(argv[1]);
        minimize_cost_function();
        //simple_grid_search();
        cleanup_globals();
    }
    catch (std::exception& exc)
    {
        cerr << endl << sep << endl
             << "Exception: " << exc.what() << endl
             << "Aborting!" << endl;
        return 1;
    }
    catch (...)
    {
        cerr << endl << sep << endl
             << "Unknown exception!" << endl
             << "Aborting!" << endl;
        return 2;
    }

    return 0;
}

