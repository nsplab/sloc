/* Test our simplex search method on various functions.
 *
 */

#include <iostream>
#include <cmath>
#include <sloc/simplex_search.h>

using namespace std;

// ----------------------------------------------------------------------------

double rosenbrock(double p[2])
{
    const double x1 = p[0];
    const double x2 = p[1];
    return 100 * pow(x2 - x1*x1, 2) + pow(1 - x1, 2);
}

double powell_quartic(double p[4])
{
    const double x1 = p[0];
    const double x2 = p[1];
    const double x3 = p[2];
    const double x4 = p[3];
    return pow(x1 + 10*x2, 2) + 5*pow(x3 - x4, 2) + pow(x2 - 2*x3, 4) + 10*pow(x1 - x4, 4);
}

double theta(double x1, double x2)
{
    double val = atan2(x2, x1);
    if (x1 < 0)
        val += M_PI;
    return val / (2 * M_PI);
}

double powell_helical_valley(double p[3])
{
    const double x1 = p[0];
    const double x2 = p[1];
    const double x3 = p[2];
    return 100 * pow(x3 - 10 * theta(x1, x2), 2) + pow(sqrt(x1*x1 + x2*x2) - 1, 2) + x3*x3;
}

// ----------------------------------------------------------------------------

void test_rosenbrock()
{

    // initial point
    double x0 = -1.2;
    double y0 = 1;

    double pt[2];
    pt[0] = x0;
    pt[1] = y0;

    sloc::NelderMead::SimplexSearch<2> simplex_search;
    simplex_search.F = rosenbrock;
    simplex_search.tol = 1e-15;
    simplex_search.debug = true;
    simplex_search.verbose = true;

    double initial_simplex[3][2];
    sloc::NelderMead::SimplexSearch<2>::make_wedge(1, pt, initial_simplex);

    double val;
    int ret = simplex_search.run(initial_simplex, pt, val);

    if (ret == 0)
    {
        cout << "Found minimum: "
             << "Point (" << pt[0] << ", " << pt[1] << ") "
             << "with value " << val << endl;
    }
}

void test_powell_quartic()
{
    // initial point
    double x0 = 3;
    double y0 = -1;
    double z0 = 0;
    double w0 = 1;

    double pt[4];
    pt[0] = x0;
    pt[1] = y0;
    pt[2] = z0;
    pt[3] = w0;

    sloc::NelderMead::SimplexSearch<4> simplex_search;
    simplex_search.F = powell_quartic;
    simplex_search.tol = 1e-15;
    simplex_search.debug = true;
    simplex_search.verbose = true;

    double initial_simplex[5][4];
    sloc::NelderMead::SimplexSearch<4>::make_wedge(1.5, pt, initial_simplex);

    double val;
    int ret = simplex_search.run(initial_simplex, pt, val);

    if (ret == 0)
    {
        cout << "Found minimum: Point ("
             << pt[0] << ", "
             << pt[1] << ", "
             << pt[2] << ", "
             << pt[3] << ") "
             << "with value " << val << endl;
    }
}

void test_powell_helical_valley()
{
    // initial point
    double x0 = -1;
    double y0 = 0;
    double z0 = 0;

    double pt[3];
    pt[0] = x0;
    pt[1] = y0;
    pt[2] = z0;

    sloc::NelderMead::SimplexSearch<3> simplex_search;
    simplex_search.F = powell_helical_valley;
    simplex_search.tol = 1e-15;
    simplex_search.debug = true;
    simplex_search.verbose = true;

    double initial_simplex[4][3];
    sloc::NelderMead::SimplexSearch<3>::make_wedge(0.8, pt, initial_simplex);

    double val;
    int ret = simplex_search.run(initial_simplex, pt, val);

    if (ret == 0)
    {
        cout << "Found minimum: Point "
             << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") "
             << "with value " << val << endl;
    }
}

// ----------------------------------------------------------------------------

int main()
{
    //test_rosenbrock();
    test_powell_quartic();
    //test_powell_helical_valley();
    return 0;
}

