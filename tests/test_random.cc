#include <iostream>
#include "random_utils.h"

int main(void)
{
    using namespace std;

    int i;
    const int N = 5;

    sloc::seed_random_generator();

    for (i = 0; i < N; i++)
        cout << "randint(2,25) = " << sloc::randint(2,25) << endl;

    for (i = 0; i < N; i++)
        cout << "uniform(0,1) = " << sloc::uniform(0,1) << endl;

    for (i = 0; i < N; i++)
        cout << "normal(5,3) = " << sloc::normal(5,3) << endl;

    return 0;
}
