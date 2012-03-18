#ifndef SLOC_NELDER_MEAD_H
#define SLOC_NELDER_MEAD_H

#include <iostream>
#include <iomanip>
#include <ctime>

namespace sloc {
namespace NelderMead {
// ----------------------------------------------------------------------------

template <int dim, int npts = dim+1>
class SimplexSearch
{
public:

    typedef double (*function_t)(double p[dim]);

    SimplexSearch()
    {
        // initialize our function to the null pointer!
        F = 0;

        // typical values for our parameters
        alpha = 1.0;
        beta = 0.5;
        gamma = 2.0;
        sigma = 0.5;

        // convergence criterion
        tol = 1e-8;

        // by default, use a high upper limit on the number of iterations
        max_iterations = 5000;

        // omit debug info by default, but print out timing information
        verbose = true;
        debug = false;
    }

    static void make_wedge(double size, double center[dim], double wedge[npts][dim])
    {
        //
        // Make a simple wedge simplex centered at point center[]
        //

        for (int i = 0; i < npts; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                if (i-1 != j)
                    wedge[i][j] = center[j] - size;
                else
                    wedge[i][j] = center[j] + size;
            }
        }
    }

    int run(double initial_simplex[npts][dim], double final_position[dim], double& final_value)
    {
        int i,j;
        double std;

        // check that the function was defined
        if (F == 0)
            return -1;

        // zero out counters
        itercount = 0;
        evalcount = 0;
        time(&t_start);
        t_tic = clock();

        if (verbose)
            std::cout << "starting simplex search" << std::endl;

        // copy the initial simplex
        for (i = 0; i < npts; i++)
        {
            for (j = 0; j < dim; j++)
            {
                P[i][j] = initial_simplex[i][j];
            }
        }


        // start by evaluating F on the simplex vertices
        for (i = 0; i < npts; i++)
        {
            y[i] = F(P[i]);
            evalcount++;

            if (verbose)
                std::cout << "function evaluation at P[" << i << "] took " << tic() << " seconds" << std::endl;
        }

        if (debug)
            print_current_simplex();

        // iterate until standard deviation of F-values at the simplex
        // vertices falls within our tolerance value
        std = std_value();
        while ((std > tol) && (itercount < max_iterations))
        {
            search_step();
            std = std_value();

            if (verbose)
                std::cout << "iteration " << itercount << " took " << tic() << " secs" << std::endl;

            if (debug)
                print_current_simplex();
        }

        // copy the results
        calc_argmaxmin();
        for (j = 0; j < dim; j++)
            final_position[j] = P[l][j];
        final_value = y[l];

        // might not have converged...
        if (itercount >= max_iterations)
        {
            return 1;
        }

        // normal return value
        return 0;
    }

private:

    void search_step()
    {
        itercount++;

        if (debug)
            std::cout << "starting iteration " << itercount << std::endl;

        calc_argmaxmin(); // sets indices (h,l) corresponding to our apex P_h and nadir P_l

        if (debug)
        {
            std::cout << "min: " << l << " "; print_point_value(P[l], y[l]); std::cout << std::endl;
            std::cout << "max: " << h << " "; print_point_value(P[h], y[h]); std::cout << std::endl;
        }

        calc_face_centroid(h); // sets P_bar using all vertices but the apex P_h

        calc_reflection(); // sets (P_star, y_star) by choosing a point away from P_h, past P_bar

        if (y_star < y[l])
        {
            //
            // the reflected value y_star is smaller than our smallest value y_l.
            //
            // let's see if we can expand again the same direction...
            // we are seeking smaller values of F, after all.
            //

            calc_expansion(); // sets (P_star2, y_star2) by expanding further away from P_bar, past P_star

            if (y_star2 < y[l])
            {
                //
                // great! the value y_star2 at our new expansion point is smaller
                // than our current nadir y_l.
                //
                // let's replace the current apex P_h with this new point P_star2.
                //
                //  P[h] = P_star2;
                //  y[h] = y_star2;
                //
                if (debug) std::cout << "P_h <- P_star2 (replace apex by expanded point)" << std::endl;
                replace_apex_by(P_star2, y_star2);
                return;
            }
            else
            {
                //
                // oops. the value y_star2 at the new expansion point P_star2 exceeded our nadir y_l.
                // let's just replace the current apex P_h with the reflected point P_star
                //
                //  P[h] = P_star;
                //  y[h] = y_star;
                //
                if (debug) std::cout << "P_h <- P_star (replace apex by reflected point)" << std::endl;
                replace_apex_by(P_star, y_star);
                return;
            }
        }
        else
        {
            //
            // ok. the value y_star at the reflected point is greater than our
            // current nadir y_l.
            //
            // so let's go down the "search gradient" along the search line
            // and try a contraction instead.
            //
            // but first, let's figure out whether the value y_star also exceeds
            // every value at non-apex vertices.
            //

            if (!every_non_apex_value_is_under(y_star))
            {
                //
                // so... there are still some non-apex vertices, say P[i],
                // whose values y[i] are bigger than y_star.
                //
                // let's be conservative and not do a contraction just yet.
                // at least not until we're closer to the apex value y_h.
                //
                // since we know that (y_l <= y_star <= y_i < y_h), we can
                // still decrease the variance of our y-values by replacing
                // our current apex P_h with its reflection P_star.
                //
                //  P[h] = P_star;
                //  y[h] = y_star;
                //
                if (debug) std::cout << "P_h <- P_star (replace apex by reflected point)" << std::endl;
                replace_apex_by(P_star, y_star);
                return;
            }
            else
            {
                //
                // ok, we are finally at the edge case where y_star is closer
                // to the apex value y_h than the values at any of the other
                // non-apex vertices. we are ready for a contraction.
                //
                // but first, let's figure out which of the two points has the
                // smallest y-value, and apply the contraction to either
                // our current simplex or its reflected cousin, depending
                // on who has the smaller apex value.
                //
                if (y_star <= y[h])
                {
                    //
                    // the reflection value y_star is lower than our current apex.
                    //
                    // let's swap P_h and P_star
                    //
                    // because we've already established that y_star was greater than all
                    // other non-apex vertices, our new P_h will still be an apex point.
                    //
                    //  P[h] = P_star;
                    //  y[h] = y_star;
                    if (debug) std::cout << "P_h <- P_star (swap apex with its reflection)" << std::endl;
                    replace_apex_by(P_star, y_star);
                }

                calc_contraction(); // sets (P_star2, y_star2) by backtracking from P_bar, back to P_h

                if (y_star2 > y[h])
                {
                    // 
                    // to summarize:
                    // - a reflection to y_star results in a value exceeding our current nadir y_l
                    // - a contraction to y_star2 results in a value exceeding our current apex y_h
                    //
                    // either way, now we aren't sure which of the two points to choose
                    //
                    // since we don't know any better, let's just move down the "search gradient"
                    // by shrinking our simplex in the direction of the point with the lowest
                    // known y-value. thus, we shift all the other vertices towards our current
                    // nadir point P_l.
                    //
                    calc_reduction(l);
                    return;
                }
                else
                {
                    //
                    // the contraction was successful! we've identified a direction in which
                    // the y-value decreases down from our current apex value y_h.
                    //
                    // since (y_l <= y_star2 <= y_h), let's decrease the variance
                    // of our y-values by replacing our current apex P_h with the
                    // contracted point P_star2.
                    //
                    if (debug) std::cout << "P_h <- P_star2 (replace apex by contraction point)" << std::endl;
                    replace_apex_by(P_star2, y_star2);
                    return;
                }
            }
        }
        return;
    }

    double total_time() const
    {
        time_t t;
        time(&t);
        return (t - t_start);
    }

    double tic()
    {
        clock_t t = clock();
        double dt = (double)(t - t_tic) / CLOCKS_PER_SEC;
        t_tic = t;
        return dt;
    }

    void print_current_simplex()
    {
        for (int i = 0; i < npts; i++)
        {
            std::cout << i << " ";
            print_point_value(P[i], y[i]);
            std::cout << std::endl;
        }
        std::cout << "time " << total_time() << std::endl;
        std::cout << "evals " << evalcount << std::endl;
        std::cout << "std " << std_value() << std::endl;
        print_sep();
    }

    void print_point_value(double pt[dim], double val)
    {
        print_point(pt);
        std::cout << " " << val;
    }

    void print_point(double pt[dim]) const
    {
        std::cout << "[ ";
        for (int j = 0; j < dim; j++)
            std::cout << pt[j] << " ";
        std::cout << "]";
    }

    void print_sep() const
    {
        for (int n = 0; n < 30; n++)
            std::cout << "-";
        std::cout << std::endl;
    }

    bool every_non_apex_value_is_under(double val)
    {
        for (int i = 0; i < npts; i++)
        {
            if (i != h) // index of apex is h
            {
                // look for a counterexample
                if (val <= y[i])
                    return false;
            }
        }
        return true;
    }

    void replace_apex_by(double pt[dim], double val)
    {
        for (int j = 0; j < dim; j++)
        {
            P[h][j] = pt[j];
        }
        y[h] = val;
    }

    double avg_value() const
    {
        double avg = 0;
        for (int i = 0; i < npts; i++)
        {
            avg += y[i];
        }
        avg /= npts;
        return avg;
    }

    double std_value() const
    {
        double ybar = avg_value();
        double yvar = 0;
        for (int i = 0; i < npts; i++)
        {
            double dy = y[i] - ybar;
            yvar += dy * dy;
        }
        yvar /= (npts - 1);
        return std::sqrt(yvar);
    }

    void calc_argmaxmin()
    {
        //
        // Calculate the respective indices (h,l) of the simplex vertices
        // corresponding to the apex point P_h and nadir point P_l,
        // relative to the function values of F at those vertices.
        //

        double yh, yl;
        h = l = 0;
        yh = yl = y[0];
        for (int i = 1; i < npts; i++)
        {
            if (y[i] > yh)
            {
                h = i;
                yh = y[i];
            }
            if (y[i] < yl)
            {
                l = i;
                yl = y[i];
            }
        }
    }

    void calc_face_centroid(int k)
    {
        //
        // Calculate the centroid of the "face" opposite vertex k.
        //

        int i,j;

        for (j = 0; j < dim; j++)
            P_bar[j] = 0;

        for (i = 0; i < npts; i++)
            for (j = 0; j < dim; j++)
                if (i != k) P_bar[j] += P[i][j];

        for (j = 0; j < dim; j++)
            P_bar[j] /= (npts-1);

        if (debug)
        {
            std::cout << "centroid: ";
            print_point(P_bar);
            std::cout << std::endl;
        }
    }

    void calc_reflection()
    {
        //
        // Choose a point P_star located away from apex P_h, aimed
        // in the direction of P_bar.
        //
        // Since P_bar is the centroid of the opposite face to the apex,
        // our new point P_star can be considered a reflection of P_h
        // across that face.
        //
        // In the linear combination below, with alpha > 0, we are
        // using an extrapolation parameter (1 + alpha) > 1.
        //

        for (int j = 0; j < dim; j++)
        {
            P_star[j] = (1 + alpha) * P_bar[j] - alpha * P[h][j];
        }

        y_star = F(P_star);
        evalcount++;

        if (debug)
        {
            std::cout << "reflection: ";
            print_point_value(P_star, y_star);
            std::cout << std::endl;
        }
    }

    void calc_expansion()
    {
        //
        // Choose a point P_star2 even further away from our apex P_h.
        //
        // We do this by expanding from the face centroid P_bar in the
        // direction of the (previously calculated) apex reflection P_star.
        //
        // In the linear combination below, gamma is an extrapolation
        // parameter. Hence gamma > 1.
        //

        for (int j = 0; j < dim; j++)
        {
            P_star2[j] = gamma * P_star[j] + (1 - gamma) * P_bar[j];
        }

        y_star2 = F(P_star2);
        evalcount++;

        if (debug)
        {
            std::cout << "expansion: ";
            std::cout << std::endl;
        }
    }

    void calc_contraction()
    {
        //
        // Choose a point P_star2 closer to our apex P_h.
        //
        // We do this by selecting a point on the line segment
        // connecting the centroid P_bar of the opposite face
        // and the original centroid.
        //
        // In the linear combination below, 0 < beta < 1, is the
        // interpolation parameter which backtracks from P_bar
        // towards the direction of the apex P_h.
        //

        for (int j = 0; j < dim; j++)
        {
            P_star2[j] = beta * P[h][j] + (1 - beta) * P_bar[j];
        }

        y_star2 = F(P_star2);
        evalcount++;

        if (debug)
        {
            std::cout << "contraction: ";
            print_point_value(P_star2, y_star2);
            std::cout << std::endl;
        }
    }

    void calc_reduction(int k)
    {
        //
        // Shrink the simplex in the direction of vertex k.
        //
        // In the reduction below we would typically use a factor
        // of sigma = 0.5, and always be shrinking in the direction
        // of the nadir point P_l.
        //

        if (debug)
        {
            std::cout << "shrinking towards " << k << " ";
            print_point(P[k]);
            std::cout << std::endl;
        }

        for (int i = 0; i < npts; i++)
        {
            if (i != k)
            {
                for (int j = 0; j < dim; j++)
                {
                    P[i][j] = P[k][j] + sigma * (P[i][j] - P[k][j]);
                }

                y[i] = F(P[i]);
                evalcount++;
            }
        }
    }

private:

    double P[npts][dim];    // current simplex
    double y[npts];         // function values at simplex vertices
    int h,l;                // indices corresponding to apex and nadir

    double P_bar[dim];      // coordinates of face centroid

    double P_star[dim];     // coordinates of reflected point
    double y_star;          // function value at reflection

    double P_star2[dim];    // coordinates of expansion/contraction point
    double y_star2;         // function value at expansion/contraction

    int itercount;
    int evalcount;
    time_t t_start;
    clock_t t_tic;

public:

    function_t F;       /// Function to minimize
    double tol;         /// convergence criterion
    double alpha;       /// reflection coefficient (alpha > 0)
    double beta;        /// contraction coefficient (0 < beta < 1)
    double gamma;       /// expansion coefficient (gamma < 1)
    double sigma;       /// reduction coefficient (0 < sigma < 1)

    int max_iterations; /// upper limit on number of iterations
    bool verbose;       /// print timing information for each iteration step
    bool debug;         /// print all debug information

};


// ----------------------------------------------------------------------------
}}
#endif
