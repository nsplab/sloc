/*
 * Parallelized version of sloc::BEM_Forward_Problem
 */

#ifndef BEM_FORWARD_PROBLEM_G_H
#define BEM_FORWARD_PROBLEM_G_H

#include <valarray>
#include <iostream>
#include <cmath>
#include "sloc/geometry.h"
#include "sloc/io_dealii.h"
#include "sloc/bem_forward_problem.h"

namespace sloc {
// ----------------------------------------------------------------------------

class BEM_Forward_Problem_G : public BEM_Forward_Problem
{
public:

    BEM_Forward_Problem_G(const Parameters& parameters);
    ~BEM_Forward_Problem_G();

    void run();
    void assemble_system();

private:
    void assemble_range_contrib(unsigned int cv_begin, unsigned int cv_end, std::valarray<double>& contrib);
    void assemble_system_from_contrib(unsigned int cv_begin, unsigned int cv_end, std::valarray<double>& contrib, const int kth_proc);
    void assemble_rhs();

public:

    // number of processes
    int numprocs;

    // current process rank
    int rank;

};

// ----------------------------------------------------------------------------
}

#endif
