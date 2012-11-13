#include <sloc/random.h>
#include <sloc/utils.h>
#include <iostream>

using namespace std;

int main(void)
{

    int i;
    const int N = 5;

    sloc::seed_random_generator();

    using sloc::randint;
    for (i = 0; i < N; i++)
        _PRINT_VALUE(randint(2,25));

    using sloc::uniform;
    for (i = 0; i < N; i++)
        _PRINT_VALUE(uniform(0,1));

    using sloc::normal;
    for (i = 0; i < N; i++)
        _PRINT_VALUE(normal(5,3));

    return 0;
}
