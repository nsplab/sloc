/*
 * (1) Read measurements Phi
 * (2) Read volume mesh.
 * (1) Calculate L.
 * (2) Calculate cost C.
 * (3) Minimize cost.
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/parameter_handler.h>
#include "bem_fwd_problem.h"
#include "dipole_sources.h"
#include "io_dealii.h"

class ForwardProblem : public sloc::BEM_ForwardProblem
{
public:
    ForwardProblem(const sloc::BEM_ForwardProblem::Parameters parameters) : sloc::BEM_ForwardProblem(parameters)
    {
    }

    void configure()
    {
        sloc::BEM_ForwardProblem::configure();
    }

    void solve(double x, double y, double z, double px, double py, double pz)
    {
        sloc::DipoleSource *src = dipole_sources._sources[0];
        src->location = dealii::Point<3>(x,y,z);
        src->dipole = dealii::Point<3>(px,py,pz);
        this->assemble_system();
        this->solve_system();
    }

    double get_phi(int m)
    {
        return phi(m);
    }
};


int main(void)
{
    using namespace std;
    using namespace dealii;

    const bool debug = true;

    const int M = 8;
    const int N = 64;

    FullMatrix<double> L;
    FullMatrix<double> Lxyz;
    FullMatrix<double> pinv;

    Vector<double> r;
    Vector<double> s;

    Vector<double> costs;
    Vector<double> phi_measured;

    Vector<double> p;
    Vector<double> p_tilde;

    int m,n;
    int i,j,k;

    int nx = 4, ny = 4, nz = 4;
    int nel = nx * ny * nz;

    double delta = 0.25;
    double x, y, z;


    L.reinit(M, N * 3);
    Lxyz.reinit(M,3);
    pinv.reinit(3,M);

    p.reinit(3);
    p_tilde.reinit(3);

    r.reinit(M);
    s.reinit(N * 3);

    costs.reinit(N);
    phi_measured.reinit(M);

    /*
     * Read phi_measured from "phi.dat"
     */
    ifstream is;
    is.open("phi.dat");
    if (!is.is_open())
    {
        cout << "Error opening 'phi.dat'" << endl;
        exit(1);
    }
    for (m = 0; m < M; m++)
        is >> phi_measured(m);
    is.close();

    /*
     * Initialize fwd problem object
     */

    deallog.depth_console(0);

    sloc::BEM_ForwardProblem::Parameters parameters;
    dealii::ParameterHandler parameter_handler;
    string prm_file = "hex.prm";

    parameters.declare_parameters(parameter_handler);
    parameter_handler.read_input(prm_file);
    parameters.get_parameters(parameter_handler);
    parameters.debug = false;
    parameters.verbose = false;

    ForwardProblem fwd(parameters);
    fwd.configure();

    /*
     * Fill in the matrix L.
     * (dimensions M x (3*N) --> 8 x 192 in this case)
     */
    for (i = 0; i < nx; i++)
    {
        x = i * delta + (delta/2);
        for (j = 0; j < ny; j++)
        {
            y = j * delta + (delta/2);
            for (k = 0; k < nz; k++)
            {
                z = k * delta + (delta/2);

                // element index from (i,j,k) triplet
                n = k * ny * nx + j * nx + i;

                // find Lx
                fwd.solve(x, y, z, 1, 0, 0);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 0) = fwd.get_phi(m);

                // find Ly
                fwd.solve(x, y, z, 0, 1, 0);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 1) = fwd.get_phi(m);

                // find Lz
                fwd.solve(x, y, z, 0, 0, 1);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 2) = fwd.get_phi(m);
            }
        }
    }

    /*
     * Calculate costs
     */
    double min_cost = numeric_limits<double>::infinity();
    int min_index = -1;
    for (n = 0; n < nel; n++)
    {
        //
        // get the L matrix for the n-th element (dimensions M x 3)
        //
        for (m = 0; m < M; m++)
        {
            Lxyz(m,0) = L(m,3*n+0);
            Lxyz(m,1) = L(m,3*n+1);
            Lxyz(m,2) = L(m,3*n+2);
        }

        //
        //              [ px ]
        // invert L --> [ py ] = pinv([Lx Ly Lz]) * phi_measured
        //              [ pz ]
        //
        pinv.left_invert(Lxyz);
        pinv.vmult(p, phi_measured);

        //
        // compute cost ||phi - L*p||
        //
        double cost2 = 0.0;
        for (m = 0; m < M; m++)
        {
            double term = phi_measured(m) - p(0) * Lxyz(m,0) - p(1) * Lxyz(m,1) - p(2) * Lxyz(m,2);
            cost2 += term * term;
        }
        costs(n) = sqrt(cost2);

        // save info associated with the best cost
        if (costs(n) < min_cost)
        {
            min_index = n;
            min_cost = costs(n);
            p_tilde = p;
        }

        if (debug)
        {
            i = n % nx;
            j = (n/nx) % ny;
            k = n / (ny * nx);
            x = i * delta + (delta/2);
            y = j * delta + (delta/2);
            z = k * delta + (delta/2);
            cout << "cost(" << n << ") = " << costs(n) << endl
                << "\tp = " << p(0) << " " << p(1) << " " << p(2) << endl
                << "\tr = " << x << " " << y << " " << z << endl
                ;
        }

    }

    if (debug)
    {
        sloc::write_vector("costs.dat", costs);
    }

    //
    // The element index is n = k * ny * nx + j * nx + i
    //
    n = min_index;
    i = n % nx;
    j = (n/nx) % ny;
    k = n / (ny * nx);
    x = i * delta + (delta/2);
    y = j * delta + (delta/2);
    z = k * delta + (delta/2);

    const string sep = "---------------------------------------------";

    // ----------------
    cout << sep << endl;
    cout << "Lxyz" << endl;
    for (m = 0; m < M; m++)
    {
        Lxyz(m,0) = L(m,3*n+0);
        Lxyz(m,1) = L(m,3*n+1);
        Lxyz(m,2) = L(m,3*n+2);
        cout << Lxyz(m,0) << " "
             << Lxyz(m,1) << " "
             << Lxyz(m,2) << endl;
    }

    // ----------------
    cout << sep << endl;
    cout << "pinv(Lxyz)" << endl;
    pinv.left_invert(Lxyz);
    for (int row = 0; row < 3; row++)
    {
        for (m = 0; m < M; m++)
            cout << pinv(row,m) << " ";
        cout << endl;
    }

    // ----------------
    cout << sep << endl;
    cout << "pinv(Lxyz) * phi_measured" << endl;
    pinv.vmult(p, phi_measured);
    cout << p(0) << endl
         << p(1) << endl
         << p(2) << endl;

    // ----------------
    cout << sep << endl;
    cout << setiosflags(ios::left);
    cout << setw(12) << "phi_measured" << "\t"
         << setw(12) << "Lxyz*p_tilde" << "\t"
         << setw(12) << "diff" << endl;
    for (m = 0; m < M; m++)
    {
        double rowprod = L(m,3*n+0) * p_tilde(0) + L(m,3*n+1) * p_tilde(1) + L(m,3*n+2) * p_tilde(2);
        double err = phi_measured(m) - rowprod;
        cout << setw(12) << phi_measured(m) << "\t"
             << setw(12) << rowprod << "\t"
             << setw(12) << err << endl;
    }

    // ----------------
    cout << sep << endl;
    cout << "Inverse problem solution" << endl;
    cout << "\tBest cost  = " << min_cost << endl;
    cout << "\t(x,y,z)    = " << x << " " << y << " " << z << endl;
    cout << "\t(px,py,pz) = " << p_tilde(0) << " " << p_tilde(1) << " " << p_tilde(2) << endl;
    cout << "\t||p||      = " << p_tilde.l2_norm() << endl;

    // ----------------

    sloc::DipoleSources actual_sources;
    actual_sources.read("data/hex.dipole");
    sloc::DipoleSource *actual_dipole_source = actual_sources._sources[0];
    Point<3> true_location = actual_dipole_source->location;
    Point<3> true_dipole = actual_dipole_source->dipole;

    cout << sep << endl;
    cout << setiosflags(ios::left);
    cout << setw(12) << "phi_measured" << "\t"
         << setw(12) << "Lxyz*p_true" << "\t"
         << setw(12) << "diff" << endl;
    double cost_with_true_p = 0.0;
    for (m = 0; m < M; m++)
    {
        double rowprod = L(m,3*n+0) * true_dipole(0) + L(m,3*n+1) * true_dipole(1) + L(m,3*n+2) * true_dipole(2);
        double err = phi_measured(m) - rowprod;
        cost_with_true_p += err * err;
        cout << setw(12) << phi_measured(m) << "\t"
             << setw(12) << rowprod << "\t"
             << setw(12) << err << endl;
    }
    cost_with_true_p = sqrt(cost_with_true_p);

    cout << sep << endl;
    cout << "phi_true" << endl;
    fwd.solve(true_location(0), true_location(1), true_location(2),
              true_dipole(0), true_dipole(1), true_dipole(2));
    for (m = 0; m < M; m++)
    {
        cout << fwd.get_phi(m) << endl;
    }

    cout << sep << endl;
    cout << "True dipole" << endl;
    cout << "\tcost       = " << cost_with_true_p << endl;
    cout << "\t(x,y,z)    = " << true_location << endl;
    cout << "\t(px,py,pz) = " << true_dipole << endl;
    cout << "\t||p||      = " << sqrt(true_dipole.square()) << endl;

    cout << sep << endl;
    cout << "Errors" << endl;
    Point<3> err_location(true_location(0) - x, true_location(1) - y, true_location(2) - z);
    Point<3> err_dipole(true_dipole(0) - p_tilde(0), true_dipole(1) - p_tilde(1), true_dipole(2) - p_tilde(2));
    double mag_err_loc = sqrt(err_location.square());
    double mag_err_dip = sqrt(err_dipole.square());
    double pct_err_loc = 100 * mag_err_loc / sqrt(true_location.square());
    double pct_err_dip = 100 * mag_err_dip / sqrt(true_dipole.square());
    cout << "\terr_location     = " << err_location << endl;
    cout << "\terr_dipole       = " << err_dipole << endl;
    cout << "\t||err_location|| = " << mag_err_loc << " (" << pct_err_loc << "%)" << endl;
    cout << "\t||err_dipole||   = " << mag_err_dip << " (" << pct_err_dip << "%)" << endl;

    return 0;
}
