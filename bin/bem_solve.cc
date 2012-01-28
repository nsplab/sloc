/* bem_solve.cc
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <deal.II/base/parameter_handler.h>
#include "bem_fwd_problem.h"

using namespace dealii;
using namespace std;

int main(int argc, char *argv[])
{
    string sep = "----------------------------------------------------";

    deallog.depth_console(3);

    deallog << "main()" << std::endl;

    try
    {
        string infile = "input.prm";
        if (argc > 1)
            infile = argv[1];

        // read parameters from file
        dealii::ParameterHandler parameter_handler;
        sloc::BEM_ForwardProblem::Parameters parameters;
        parameters.declare_parameters(parameter_handler);
        parameter_handler.read_input(infile);
        parameter_handler.print_parameters(std::cout, dealii::ParameterHandler::ShortText);
        parameters.get_parameters(parameter_handler);

        // initialize forward problem and run it
        sloc::BEM_ForwardProblem fwd_problem(parameters);
        fwd_problem.run();
    }
    catch (std::exception &exc)
    {
        cerr << endl
             << sep << endl
             << "Exception on processing: " << endl
             << exc.what() << endl
             << "Aborting!" << endl
             << sep << endl;
        return 1;
    }
    catch (...)
    {
        cerr << endl
             << sep << endl
             << "Unknown exception!" << endl
             << "Aborting!" << endl
             << sep << endl;
        return 1;
    }

    return 0;
}
