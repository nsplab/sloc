/* bem_solve.cc
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include "bem_fwd_problem.h"

using namespace dealii;
using namespace std;

int main(void)
{
    string sep = "----------------------------------------------------";

    deallog.depth_console(3);

    deallog << "main()" << std::endl;

    try
    {
        sloc::BEM_ForwardProblem fwd_problem;
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
