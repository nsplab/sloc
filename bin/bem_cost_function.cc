/* bem_cost_function.cc
 *
 * Evaluate forward solution's cost function in a cloud of points around the true solution.
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <valarray>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/smart_ptr.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include <sloc/bem_forward_problem.h>
#include <sloc/bem_forward_problem_par.h>
#include <sloc/io_dealii.h>

namespace fs = boost::filesystem;

using namespace std;
using boost::shared_ptr;

// ----------------------------------------------------------------------------

class ForwardProblem : public sloc::BEM_Forward_Problem_P
{
public:

    ForwardProblem(const sloc::BEM_Forward_Problem::Parameters& fwd_params)
        : sloc::BEM_Forward_Problem_P(fwd_params)
    {
    }

    void configure()
    {
        // first, create a dummy dipole source
        dealii::Point<3> dummy_x(0,0,0), dummy_p(0,0,1);
        dipole_sources.add_source(dummy_x, dummy_p);

        // configure the rest of the parameters
        sloc::BEM_Forward_Problem_P::configure();
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

class CostFunctionParameters
{
public:
    static void declare_parameters(dealii::ParameterHandler& prm)
    {
        using namespace dealii;

        prm.declare_entry("electrodes", "electrodes.dat", Patterns::Anything(), "Input file with electrode potential measurements");
        prm.declare_entry("output_file", "cost.dat", Patterns::Anything(), "Values of cost function on all grid points");
        //prm.declare_entry("output_sources", "sources.dipole", Patterns::Anything(), "Final output of source localization");
        //prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level of inverse solver");
        //prm.declare_entry("debug", "false", Patterns::Bool(), "Debugging information for inverse solver");

        prm.enter_subsection("Forward Problem Parameters");
        {
            prm.declare_entry("surface_mesh", "", Patterns::Anything(), "Filename for surface mesh");
            prm.declare_entry("surface_mesh_materials", "", Patterns::Anything(), "Filename for surface mesh materials");
            prm.declare_entry("material_data", "", Patterns::Anything(), "Filename for material data");
            prm.declare_entry("verbose", "true", Patterns::Bool(), "Verbosity level of forward solver");
            prm.declare_entry("debug", "false", Patterns::Bool(), "Debugging information for forward solver");
        }
        prm.leave_subsection();

        prm.enter_subsection("Grid Parameters");
        {
            prm.declare_entry("grid_dims", "10,10,10", Patterns::List(Patterns::Integer(), 3, 3), "Number of subdivisions along each axis of the grid");
            prm.declare_entry("grid_lengths", "0.01,0.01,0.01", Patterns::List(Patterns::Double(), 3, 3), "Length of grid along each axis in meters");
            prm.declare_entry("grid_center", "0,0,0", Patterns::List(Patterns::Double(), 3, 3), "Center of grid");
        }
        prm.leave_subsection();
    }

    void get_parameters(dealii::ParameterHandler& prm)
    {
        electrodes = prm.get("electrodes");
        output_file = prm.get("output_file");
        //output_sources = prm.get("output_sources");
        //verbose = prm.get_bool("verbose");
        //debug = prm.get_bool("debug");

        prm.enter_subsection("Forward Problem Parameters");
        {
            surface_mesh = prm.get("surface_mesh");
            surface_mesh_materials = prm.get("surface_mesh_materials");
            material_data = prm.get("material_data");
            fwd_verbose = prm.get_bool("verbose");
            fwd_debug = prm.get_bool("debug");
        }
        prm.leave_subsection();

        prm.enter_subsection("Grid Parameters");
        {
            using namespace boost::algorithm;
            string str;

            str = prm.get("grid_dims");
            vector<string> dims;
            split(dims, str, is_any_of(", "), token_compress_on);
            grid_dims[0] = boost::lexical_cast<int>(dims[0]);
            grid_dims[1] = boost::lexical_cast<int>(dims[1]);
            grid_dims[2] = boost::lexical_cast<int>(dims[2]);

            str = prm.get("grid_lengths");
            vector<string> lengths;
            split(lengths, str, is_any_of(", "), token_compress_on);
            grid_lengths[0] = boost::lexical_cast<double>(lengths[0]);
            grid_lengths[1] = boost::lexical_cast<double>(lengths[1]);
            grid_lengths[2] = boost::lexical_cast<double>(lengths[2]);

            str = prm.get("grid_center");
            vector<string> center;
            split(center, str, is_any_of(", "), token_compress_on);
            grid_center[0] = boost::lexical_cast<double>(center[0]);
            grid_center[1] = boost::lexical_cast<double>(center[1]);
            grid_center[2] = boost::lexical_cast<double>(center[2]);
        }
        prm.leave_subsection();

    }

    void print() const
    {
        using namespace std;
        cout << "electrodes = " << electrodes << endl;
        cout << "output_file = " << output_file << endl;
        cout << "surface_mesh = " << surface_mesh << endl;
        cout << "surface_mesh_materials = " << surface_mesh_materials << endl;
        cout << "material_data = " << material_data << endl;
        cout << "grid_dims = " << grid_dims[0] << ", " << grid_dims[1] << ", " << grid_dims[2] << endl;
        cout << "grid_lengths = " << grid_lengths[0] << ", " << grid_lengths[1] << ", " << grid_lengths[2] << endl;
        cout << "grid_center = " << grid_center[0] << ", " << grid_center[1] << ", " << grid_center[2] << endl;
        //cout << "verbose = " << verbose << endl;
        //cout << "debug = " << debug << endl;
        cout << "fwd_verbose = " << fwd_verbose << endl;
        cout << "fwd_debug = " << fwd_debug << endl;
    }

public:
    std::string electrodes;
    std::string output_file;
    std::string surface_mesh;
    std::string surface_mesh_materials;
    std::string material_data;
    int grid_dims[3];
    double grid_lengths[3];
    double grid_center[3];
    //bool verbose, debug;
    bool fwd_verbose, fwd_debug;
};

class CostFunction
{
private:
    shared_ptr<ForwardProblem> forward_problem;
    sloc::BEM_Forward_Problem::Parameters fwd_prm;

public:

    CostFunction(CostFunctionParameters& parameters)
    {
        // Initialize forward problem
        fwd_prm.surface_mesh = parameters.surface_mesh;
        fwd_prm.surface_mesh_materials = parameters.surface_mesh_materials;
        fwd_prm.material_data = parameters.material_data;
        fwd_prm.dipole_sources = "";
        fwd_prm.verbose = parameters.fwd_verbose;
        fwd_prm.debug = parameters.fwd_debug;
        fwd_prm.logfile = "";

        // Initialize cost function
        forward_problem = shared_ptr<ForwardProblem>(new ForwardProblem(fwd_prm));

        // disable dealii logging (fwd solver needs to be quiet)
        dealii::deallog.depth_console(0);

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
        forward_problem->configure();

        int m;
        const int M = num_electrodes;

        double x = pt[0];
        double y = pt[1];
        double z = pt[2];

        // find L[:,0]
        forward_problem->solve(x, y, z, 1, 0, 0);
        for (m = 0; m < M; m++)
            L(m,0) = forward_problem->get_phi(electrode_dof_index[m]);

        // find L[:,1]
        forward_problem->solve(x, y, z, 0, 1, 0);
        for (m = 0; m < M; m++)
            L(m,1) = forward_problem->get_phi(electrode_dof_index[m]);

        // find L[:,2]
        forward_problem->solve(x, y, z, 0, 0, 1);
        for (m = 0; m < M; m++)
            L(m,2) = forward_problem->get_phi(electrode_dof_index[m]);

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

public:

    unsigned int num_electrodes;
    std::vector<unsigned int> electrode_dof_index;
    dealii::Vector<double> phi_measured;

    dealii::FullMatrix<double> L;
    dealii::FullMatrix<double> L_pinv;
    dealii::Vector<double> dipole;

};
// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();

    string sep = "----------------------------------------------------";

    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <parameter-file>" << endl;
        return 0;
    }

    if (!fs::exists(argv[1]))
    {
        cerr << "Could not find parameter file '" << argv[1] << "'" << endl;
        return 1;
    }
    string parameter_file = argv[1];

    if (rank == 0)
    {
        if (argc < 2)
        {
            cout << "Usage: " << argv[0] << " <input.prm>" << endl;
            MPI::Finalize();
            return 1;
        }

    }

    try
    {
        CostFunctionParameters parameters;
        dealii::ParameterHandler parameter_handler;
        parameters.declare_parameters(parameter_handler);
        parameter_handler.read_input(parameter_file);
        parameters.get_parameters(parameter_handler);

        if (rank == 0)
        {
            parameters.print();
        }

        CostFunction cost_function(parameters);

        //
        // Build grid and evaluate cost function
        //
        int i, j, k;

        int nx = parameters.grid_dims[0];
        int ny = parameters.grid_dims[1];
        int nz = parameters.grid_dims[2];

        double dx = parameters.grid_lengths[0] / nx;
        double dy = parameters.grid_lengths[1] / ny;
        double dz = parameters.grid_lengths[2] / nz;

        double center[3];
        center[0] = parameters.grid_center[0];
        center[1] = parameters.grid_center[1];
        center[2] = parameters.grid_center[2];

        double x0, y0, z0;
        x0 = center[0] - (nx * dx / 2);
        y0 = center[1] - (ny * dy / 2);
        z0 = center[2] - (nz * dz / 2);

        valarray<double> X(nx);
        for (i = 0; i < nx; i++)
            X[i] = x0 + i * dx;

        valarray<double> Y(ny);
        for (j = 0; j < ny; j++)
            Y[j] = y0 + j * dy;

        valarray<double> Z(nz);
        for (k = 0; k < nz; k++)
            Z[k] = z0 + k * dz;

        ofstream out(parameters.output_file.c_str());
        for (i = 0; i < nx; i++)
        {
            double x = X[i];
            for (j = 0; j < ny; j++)
            {
                double y = Y[j];
                for (k = 0; k < nz; k++)
                {
                    double z = Z[k];

                    double pt[3] = {x,y,z};
                    double dp[3] = {0,0,0};
                    double cost = cost_function.cost_for_single_source_at(pt, dp);

                    out << pt[0] << " " << pt[1] << " " << pt[2] << " "
                        << dp[0] << " " << dp[1] << " " << dp[2] << " "
                        << cost << endl;
                }
            }
        }
        out.close();
        
    }
    catch (std::exception& exc)
    {
        cerr << endl << sep << endl
             << "Exception: " << exc.what() << endl
             << "Aborting!" << endl;
        MPI::Finalize();
        return 1;
    }
    catch (...)
    {
        cerr << endl << sep << endl
             << "Unknown exception!" << endl
             << "Aborting!" << endl;
        MPI::Finalize();
        return 1;
    }

    MPI::Finalize();

    return 0;
}
