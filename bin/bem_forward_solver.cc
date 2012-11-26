// bem_forward_solver.cc

#include <iostream>
#include <fstream>
#include <string>

#include <boost/filesystem.hpp>
#include <deal.II/base/parameter_handler.h>
#include <sloc/bem_forward_problem.h>
#include <sloc/bem_forward_problem_par.h>

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace dealii;
    namespace fs = boost::filesystem;

    string sep = "----------------------------------------------------";

    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    int numprocs = MPI::COMM_WORLD.Get_size();

    if (rank == 0)
    {
        if (argc < 2)
        {
            cout << "Usage: " << argv[0] << " <input.prm>" << endl;
            MPI::Finalize();
            return 1;
        }

    }

    if (!fs::exists(argv[1]))
    {
        cout << "Error: parameter file '" << argv[1] << "' does not exist!" << endl;
        MPI::Finalize();
        return 2;
    }

    string infile = argv[1];

    try
    {
        // read parameters from file
        dealii::ParameterHandler parameter_handler;
        sloc::BEM_Forward_Problem::Parameters parameters;
        parameters.declare_parameters(parameter_handler);
        parameter_handler.read_input(infile);
        parameters.get_parameters(parameter_handler);

        // logging settings
        if (parameters.debug && (rank == 0))
        {
            deallog.depth_console(3);
            deallog << "main()" << endl;
            parameter_handler.print_parameters(std::cout, dealii::ParameterHandler::ShortText);
        }
        else
        {
            deallog.depth_console(0);
        }

        // silence all nonzero processes, regardless of what the .prm file says
        if (rank > 0)
        {
            parameters.debug = false;
            parameters.verbose = false;
        }

        // initialize forward problem and run it
        if (numprocs > 1)
        {
            sloc::BEM_Forward_Problem_P fwd_problem(parameters);
            fwd_problem.run();
        }
        else
        {
            sloc::BEM_Forward_Problem fwd_problem(parameters);
            fwd_problem.run();
        }
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

