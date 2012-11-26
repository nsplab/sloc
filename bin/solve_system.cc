/*
 * Read "system_matrix.dat" and "system_rhs.dat" and solve the linear system they represent.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <gmm/gmm.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <sloc/io_dealii.h>
#include <sloc/utils.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

// ----------------------------------------------------------------------------

void read_std_vector(const char *filename, std::vector<double>& vec)
{
    ifstream in;
    in.open(filename);

    if (!in)
    {
        std::stringstream ss;
        ss << "Invalid filename: " << filename;
        throw std::runtime_error(ss.str());
    }

    vec.clear();
    while (!in.eof())
    {
        double num;
        in >> num;
        if (!in.eof())
            vec.push_back(num);
    }

    in.close();
}

void read_vector(const char *filename, dealii::Vector<double>& vec)
{
    std::vector<double> v;
    read_std_vector(filename, v);

    const int n = gmm::vect_size(v);
    vec.reinit(n);

    // copy into dealii object
    for (int i = 0; i < n; ++i)
        vec(i) = v[i];
}

void read_matrix(const char *filename, dealii::FullMatrix<double>& mat)
{
    std::vector<double> v;
    read_std_vector(filename, v);

    // XXX: assert that v_size is the square of an integer
    const int n = static_cast<int>(floor(sqrt(gmm::vect_size(v))));
    mat.reinit(n,n);

    // copy into dealii object
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat(i,j) = v[j + n * i];
}

// ----------------------------------------------------------------------------

void run_dealii_linear_solver(const po::variables_map& vm)
{
    using namespace dealii;
    stringstream ss;

    // command line arguments
    string mfile = vm["matrix"].as<string>();
    string rfile = vm["rhs"].as<string>();
    string outfile = vm["output"].as<string>();

    // represents matrix of linear system: Ax = b
    dealii::FullMatrix<double> system_matrix;
    read_matrix(mfile.c_str(), system_matrix);
    if (system_matrix.m() != system_matrix.n())
    {
        ss << "Matrix of linear system is not square! "
           << "Dimensions are (" << system_matrix.m() << " x " << system_matrix.n() << ")";
        throw std::invalid_argument(ss.str());
    }
    cout << "Read (" << system_matrix.m() << " x " << system_matrix.n() << ") matrix from "
         << fs::path(mfile) << endl;

    // represents right hand side of linear system: Ax = b
    dealii::Vector<double> system_rhs;
    read_vector(rfile.c_str(), system_rhs);
    if (system_matrix.m() != system_rhs.size())
    {
        ss << "Dimensions of matrix and RHS vector do not match!";
        throw std::invalid_argument(ss.str());
    }
    cout << "Read (" << system_rhs.size() << " x 1) column vector from "
         << fs::path(rfile) << endl;

    const unsigned int n = system_rhs.size();

    // represents solution of linear system
    dealii::Vector<double> x;
    x.reinit(n);

    // measures wallclock time
    Timer timer;

    // controls linear solver (tolerance, logging, etc.)
    SolverControl solver_control;
    solver_control.log_frequency(vm["log-frequency"].as<int>());
    solver_control.log_history(vm["log-history"].as<bool>());
    solver_control.log_result(vm["log-result"].as<bool>());
    solver_control.set_max_steps(vm["max-steps"].as<int>());
    solver_control.set_tolerance(vm["tolerance"].as<double>());

    // solve the linear system!
    SolverGMRES< Vector<double> > solver(solver_control);
    solver.solve(system_matrix, x, system_rhs, dealii::PreconditionIdentity());

    // save the solution
    sloc::write_vector(outfile.c_str(), x);
    cout << "Wrote solution to " << fs::path(outfile) << endl;

    // done!
}

void run_getfem_linear_solver(const po::variables_map& vm)
{
    throw std::runtime_error("procedure not yet implemented!");
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    try
    {
        // process command line args
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "produce help message")
            ("output,o", po::value<string>(), "output file for solution to linear system")
            ("matrix,m", po::value<string>()->default_value("system_matrix.dat"), "matrix of linear system")
            ("rhs,r", po::value<string>()->default_value("system_rhs.dat"), "rhs of linear system")
            ("log-frequency", po::value<int>()->default_value(1), "log frequency for solver")
            ("log-history", po::value<bool>()->default_value(false)->zero_tokens(), "log history in solver")
            ("log-result", po::value<bool>()->default_value(true)->zero_tokens(), "log result in solver")
            ("max-steps", po::value<int>()->default_value(100), "max steps for iterative solver")
            ("tolerance", po::value<double>()->default_value(1e-10), "convergence tolerance")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help") || (argc == 1))
        {
            cout << argv[0] << " OPTIONS\n\n";
            cout << desc << endl;
            return 0;
        }

        if (!vm.count("output"))
        {
            cout << "Required flag '--output' is missing!" << endl;
            return 1;
        }

        if (!fs::exists(vm["matrix"].as<string>()))
        {
            cerr << "Matrix input file " << vm["matrix"].as<string>() << " does not exist!" << endl;
            return 1;
        }

        if (!fs::exists(vm["rhs"].as<string>()))
        {
            cerr << "RHS input file " << vm["rhs"].as<string>() << " does not exist!" << endl;
            return 1;
        }

        // run our procedure!
        run_dealii_linear_solver(vm);
    }
    catch (std::exception& exc)
    {
        cerr << "Exception: " << exc.what() << endl;
        cerr << "Aborting!" << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

    return 0;
}

