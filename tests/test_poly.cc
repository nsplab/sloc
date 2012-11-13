// test_poly.cc


#include <sloc/utils.h>
#include <gmm/gmm_except.h>

// Basic Geometric Tools
#include <getfem/bgeot_config.h>
#include <getfem/bgeot_poly.h>

using namespace std;
using bgeot::short_type;

void test_poly()
{
    bgeot::base_poly W, Z;

    W[0] = 1.0;
    Z[0] = 2.0;

    _PRINT_VALUE(W.real_degree());
    _PRINT_VALUE(W);
    _PRINT_EXEC(W.direct_product(Z));
    _PRINT_VALUE(W);
}

int main(int argc, const char *argv[])
{
    GMM_SET_EXCEPTION_DEBUG;    // Exceptions make a memory fault, to debug
    FE_ENABLE_EXCEPT;           // Enable floating point exception for NaN

    try {
        test_poly();
    } GMM_STANDARD_CATCH_ERROR;

    return 0;
}
