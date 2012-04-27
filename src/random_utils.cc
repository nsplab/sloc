#include "random_utils.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>

static boost::mt19937 gen;

void sloc::seed_random_generator()
{
    gen.seed(static_cast<unsigned int>(time(0)));
}

int sloc::randint(int a, int b)
{
    boost::uniform_int<> dist(a,b);
    return dist(gen);
}

double sloc::uniform(double a, double b)
{
    boost::uniform_int<> dist(a,b);
    return dist(gen);
}

double sloc::normal(double mu, double sigma)
{
    boost::normal_distribution<> dist(mu, sigma);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > N(gen, dist);
    return N();
}

// EOF
