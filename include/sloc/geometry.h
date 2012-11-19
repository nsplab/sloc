#ifndef SLOC_GEOMETRY_H
#define SLOC_GEOMETRY_H

#include "sloc/random.h"
#include <deal.II/base/point.h>
#include <getfem/bgeot_vector.h>

namespace sloc
{
    template <class PT>
    void generate_random_points(int N, std::vector<PT>& xs, PT& pmin, PT& pmax)
    {
        double x,y,z;
        double mu = 0;
        double sigma = 1;
        for (int i = 0; i < N; i++)
        {
            x = sloc::normal(mu,sigma);
            y = sloc::normal(mu,sigma);
            z = sloc::normal(mu,sigma);
            PT p(x,y,z);
            xs.push_back(p);
            //bbox.update(p);
        }
    }

    dealii::Point<3> cross(dealii::Point<3> A, dealii::Point<3> B);

    dealii::Point<3> triangle_normal(const bgeot::base_matrix& G);

}

#endif
