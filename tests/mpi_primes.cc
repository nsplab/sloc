// mpi_primes.cc
// Based on http://people.sc.fsu.edu/~jburkardt/cpp_src/prime_mpi/prime_mpi.html
//
#include "config.h"
#include <iostream>
#include <iomanip>

#if defined(USING_MPI)
// ----------------------------------------------------------------------------
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "mpi.h"

using namespace std;

void timestamp()
{
    const int TIME_SIZE = 40;
    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);
    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    cout << time_buffer << endl;
    return;
}

int prime_number(int n, int id, int p)
{
    //
    // Outputs the number of prime numbers up to N
    //
    // In order to divide the work evenly among P processors,
    // processor ID starts at 2+ID and skips by P.
    //
    // A naive algorithm is used.
    //
    // Input N      the maximum number to check
    // Input ID     the ID of this process, between 0 and P-1
    // Input P      the number of processes
    //
    int i,j;
    int prime;
    int total;

    total = 0;
    for (i = 2 + id; i <= n; i += p)
    {
        prime = 1;
        for (j = 2; j < i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

int main(int argc, char *argv[])
{
    int id;
    int n, n_factor, n_hi, n_lo;
    int p, primes, primes_part;
    double wtime;

    n_lo = 1;
    n_hi = 131072;
    n_factor = 2;

    // Initialize MPI
    MPI::Init(argc, argv);

    // Get the number of processes
    p = MPI::COMM_WORLD.Get_size();

    // Determine this processes' rank:w
    id = MPI::COMM_WORLD.Get_rank();

    if (id == 0)
    {
        timestamp();
        cout << "\n"
             << "MPI_PRIMES\n\n"
             << "  An MPI program to count the number of primes.\n"
             << "  The number of processes is " << p << "\n\n"
             << "         N        Pi          Time\n"
             << endl;
    }

    for (n = n_lo; n <= n_hi; n *= n_factor)
    {
        if (id == 0)
        {
            wtime = MPI::Wtime();
        }

        MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0);

        primes_part = prime_number(n, id, p);

        MPI::COMM_WORLD.Reduce(&primes_part, &primes, 1, MPI::INT, MPI::SUM, 0);

        if (id == 0)
        {
            wtime = MPI::Wtime() - wtime;

            cout << "  " << setw(8) << n
                 << "  " << setw(8) << primes
                 << "  " << setw(14) << wtime
                 << endl;
        }
    }

    // Terminate MPI
    MPI::Finalize();

    // Terminate
    if (id == 0)
    {
        cout << "\n"
             << "MPI_PRIMES - Master process:\n"
             << "  Normal end of execution.\n"
             << endl;
        timestamp();
    }

    return 0;
}

#else
// ----------------------------------------------------------------------------

int main(void)
{
    using namespace std;
    cerr << "MPI not available" << endl;
    return 1;
}

#endif // USING_MPI
