#include "sloc/geometry.h"

dealii::Point<3> sloc::cross(dealii::Point<3> A, dealii::Point<3> B)
{
    //
    //             | i    j   k |
    // C = A x B = | a0  a1  a2 |
    //             | b0  b1  b2 |
    //
    //   = (a1*b2 - a2*b1)i - (a0*b2 - a2*b0)j + (a0*b1 - a1*b0)k
    //
    double C[3];
    C[0] = +(A(1)*B(2) - A(2)*B(1));
    C[1] = -(A(0)*B(2) - A(2)*B(0));
    C[2] = +(A(0)*B(1) - A(1)*B(0));
    return dealii::Point<3>(C[0], C[1], C[2]);
} 

dealii::Point<3> sloc::triangle_normal(const bgeot::base_matrix& G)
{
    // columns of G are points in the corresponding simplex
    dealii::Point<3> A(G(0,0), G(1,0), G(2,0));
    dealii::Point<3> B(G(0,1), G(1,1), G(2,1));
    dealii::Point<3> C(G(0,2), G(1,2), G(2,2));
    dealii::Point<3> N = cross(C-A, B-A);
    return N / std::sqrt(N.square());
}

