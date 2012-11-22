/*
 * Parallelized version of sloc::BEM_Forward_Problem
 */

#ifndef BEM_FORWARD_PROBLEM_PAR_H
#define BEM_FORWARD_PROBLEM_PAR_H

#include <valarray>
#include <iostream>
#include <cmath>
#include "sloc/geometry.h"
#include "sloc/io_dealii.h"
#include "sloc/bem_forward_problem.h"

namespace sloc {
// ----------------------------------------------------------------------------

class BEM_Forward_Problem_P : public BEM_Forward_Problem
{
public:

    BEM_Forward_Problem_P(const Parameters& parameters);
    ~BEM_Forward_Problem_P();

    void run();
    void assemble_system();

public:

    // number of processes
    int numprocs;

    // current process rank
    int rank;

};

// ----------------------------------------------------------------------------
}

#endif