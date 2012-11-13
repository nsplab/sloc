#ifndef SLOC_RANDOM_H
#define SLOC_RANDOM_H

namespace sloc
{
    void seed_random_generator();
    int randint(int a, int b);
    double uniform(double a, double b);
    double normal(double mu, double sigma);
}

#endif
