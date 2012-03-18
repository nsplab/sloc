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
#include <deal.II/lac/full_matrix.h>
#include "bem_fwd_problem.h"
#include "dipole_sources.h"
#include "io_dealii.h"
#include "nelder_mead.h"

// ----------------------------------------------------------------------------

class ForwardProblem : public sloc::BEM_ForwardProblem
{

public:

    ForwardProblem(const sloc::BEM_ForwardProblem::Parameters parameters)
        : sloc::BEM_ForwardProblem(parameters)
    {
    }

    void configure()
    {
        sloc::BEM_ForwardProblem::configure();
    }

    void solve(double x, double y, double z, double px, double py, double pz)
    {
        // Assume we only have a single dipole
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

    InverseProblem(const char *parameters_file)
    {
        //
        // Initialize parameters object for forward-problem
        //

        dealii::deallog.depth_console(0);

        dealii::ParameterHandler parameter_handler;
        parameters.declare_parameters(parameter_handler);

        std::string filename = parameters_file;
        parameter_handler.read_input(filename);
        parameters.get_parameters(parameter_handler);

        //
        // Initialize electrodes
        //

        Assert(!parameters.electrodes.empty(), dealii::ExcEmptyObject());
        std::ifstream is;
        is.open(parameters.electrodes.c_str());

        is >> num_electrodes;
        const int M = num_electrodes;

        unsigned int code;
        is >> code;

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

    double cost_for_single_source_at(double pt[3])
    {
        ForwardProblem forward_problem(parameters);
        forward_problem.configure();

        const int M = num_electrodes;
        int m;

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

    double cost_for_sources_at(double *positions)
    {
        //
        // XXX: implement this method
        //
        double cost = 0.0;
        return cost ;
    }


private:

    sloc::BEM_ForwardProblem::Parameters parameters;

    unsigned int num_electrodes;
    std::vector<unsigned int> electrode_dof_index;
    dealii::Vector<double> phi_measured;

    dealii::FullMatrix<double> L;
    dealii::FullMatrix<double> L_pinv;
    dealii::Vector<double> dipole;

};

/* global object! */
InverseProblem *g_inverse_problem = 0;

// ----------------------------------------------------------------------------

void init_inverse_problem(char *parameters_file)
{
    g_inverse_problem = new InverseProblem(parameters_file);
}

double cost_function(double pt[3])
{
    return g_inverse_problem->cost_for_single_source_at(pt);
}

void minimize_cost_function()
{
    const int dim = 3;
    double initial_point[dim] = {0.5, 0.5, 0.5};
    double neighborhood = 0.2;

    sloc::NelderMead::SimplexSearch<dim> simplex_search;
    simplex_search.F = cost_function;
    simplex_search.tol = 1e-8;
    simplex_search.debug = true;
    simplex_search.verbose = true;

    double initial_simplex[4][dim];
    sloc::NelderMead::SimplexSearch<dim>::make_wedge(neighborhood, initial_point, initial_simplex);

    double final_cost;
    double final_point[dim];
    int ret = simplex_search.run(initial_simplex, final_point, final_cost);

    if (ret == 0)
    {
        std::cout << "Found minimum at ("
                  << final_point[0] << ", "
                  << final_point[1] << ", "
                  << final_point[2] << ") with cost "
                  << final_cost
                  << std::endl;
    }

    // XXX: write out dipoles file
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    if (argc > 1)
    {
        init_inverse_problem(argv[1]);
        minimize_cost_function();
    }
    else
    {
        std::cout << "Usage: " << argv[0] << " input.prm" << std::endl;
        return -1;
    }

    return 0;
}

