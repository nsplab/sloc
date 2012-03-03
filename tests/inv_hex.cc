/* 
 * (1) Read measurements Phi
 * (2) Read volume mesh.
 * (1) Calculate L.
 * (2) Calculate cost C.
 * (3) Minimize cost.
 *
 */
#include <iostream>
#include <string>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include "bem_fwd_problem.h"


struct ForwardProblem
{
public:
    void solve(double x, double y, double z, double px, double py, double pz)
    {
    }
public:
    dealii::Vector<int> dofs;
    dealii::Vector<double> phi;
};

void find_optimal_p(
    const dealii::Vector<double>& phi,
    const dealii::Vector<double>& Lx,
    const dealii::Vector<double>& Ly,
    const dealii::Vector<double>& Lz,
    double& px,
    double& py,
    double& pz
    )
{
    // (1) calculate pseudo inverse of [Lx Ly Lz]: pinv(Le) = inv(Le' * Le) * Le'
    // (2) multiply pinv(L) by phi
}

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace dealii;

    FullMatrix<double> L;
    Vector<double> s;
    Vector<double> r;

    Vector<double> Lx, Ly, Lz;
    Vector<double> costs;
    Vector<double> phi_measured;
    
    double meas[8] = {
        0.589543,
        -1.16486,
        -1.69593,
        -1.8435,
        27.126,
        -4.5942,
        -2.14294,
        -3.57749
    };

    const int M = 8;
    const int N = 64;


    phi_measured.reinit(M);
    for (m = 0; m < M; m++)
        phi_measured(m) = meas[m];

    r.reinit(M);
    s.reinit(N * 3);
    L.reinit(M, N * 3);

    Lx.reinit(M);
    Ly.reinit(M);
    Lz.reinit(M);
    costs.reinit(N);

    int m,n;
    int i,j,k;
    int a,b,c;


    int nx = 4, ny = 4, nz = 4;
    double delta = 0.25;
    double x, y, z;

    double px, py, pz;


    /*
     * Fill in the matrix L. (dimensions 8 x 192)
     */
    n = 0;
    for (i = 0; i < nx; i++)
    {
        x = i * delta + (delta/2);
        for (j = 0; j < ny; j++)
        {
            y = j * delta + (delta/2);
            for (k = 0; k < nz; k++)
            {
                z = k * delta + (delta/2);

                // find Lx
                fwd.solve(x, y, z, 1, 0, 0);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 0) = fwd.phi(m)

                // find Ly
                fwd.solve(x, y, z, 0, 1, 0);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 1) = fwd.phi(m);

                // find Lz
                fwd.solve(x, y, z, 0, 0, 1);
                for (m = 0; m < M; m++)
                    L(m, 3*n + 2) = fwd.phi(m);

                // go to next element index
                n++;
            }
        }
    }

    /* 
     * Calculate cost
     */
    n = 0;
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            for (k = 0; k < nz; k++)
            {
                // get L columns for our element n
                for (m = 0; m < M; m++)
                {
                    Lx(m) = L(m,3*n+0);
                    Ly(m) = L(m,3*n+1);
                    Lz(m) = L(m,3*n+2);
                }

                //
                // invert L --> [px py pz]' = pinv([Lx Ly Lz]) * phi_measured
                //
                find_optimal_p(phi, Lx, Ly, Lz, px, py, pz);

                // compute cost
                cost(n) = 0.0;
                for (m = 0; m < M; m++)
                {
                    double term = phi(m) - px * Lx(m) - py * Ly(m) - pz * Lz(m);
                    cost(n) += term * term;
                }
                cost(n) = sqrt(cost(n));

                // next element index
                n++;
            }
        }
    }


    return 0;
}
