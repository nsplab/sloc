/*
 * Test MPI_Send and MPI_Recv on array data
 */
#include <sloc/config.h>
#include <sloc/utils.h>
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{ 
    // initialize MPI
    MPI::Init(argc, argv);

    // get the number of processes
    int P = MPI::COMM_WORLD.Get_size();

    // determine this process' rank
    int p = MPI::COMM_WORLD.Get_rank();

    // size of arrays
    int N = 1000;
    int n = N / P;

    // the array and subarray
    double *A = 0;
    double *Asub = 0;
    int i;

    // first, print some info
    if (p == 0)
    {
        _PRINT_VALUE(P);
        _PRINT_VALUE(N);
        _PRINT_VALUE(N / P);
        _PRINT_VALUE(N % P);
    }

    /*
     * Only process 0 needs to know about the full array A. Everyone else
     * works with the subarray Asub.
     */

    if (p == 0)
    {
        A = new double[N];
    }

    Asub = new double[n];

    /*
     * Now, let all processes work on their respective subarray.
     * Note that later, you must use a stride of n when copying back into A.
     */
    for (i = 0; i < n; i++)
    {
        Asub[i] = p;
    }

    /*
     * Send back results to process 0. 
     *
     * Note that process 0 itself handles the indices in the intervals
     * [0,n) and [n*P,N), whereas the processes k = 1,...,P-1 are responsible
     * for the indices in the interval [k*n,(k+1)*n)
     */
    int dest, source, tag;
    MPI::Status status;

    if (p != 0)
    {
        // send subarray to process 0
        dest = 0; tag = 0;
        MPI::COMM_WORLD.Send(Asub, n, MPI::DOUBLE, dest, tag);
    }
    else
    {
        // copy own data to A directly (no need to receive from self)
        for (i = 0; i < n; i++)
            A[i] = Asub[i];

        // take care of remainder (again, directly into array A)
        for (i = n*P; i < N; i++)
            A[i] = P;

        // receive everything from the other processes
        for (int k = 1; k < P; k++)
        {
            source = k; tag = 0;
            MPI::COMM_WORLD.Recv(Asub, n, MPI::DOUBLE, source, tag, status);

            // copy subarray Asub into A
            for (i = 0; i < n; i++)
                A[i + k * n] = Asub[i];
        }

    }

    /*
     * Finally, let process 0 do its thing with array A.
     */
    if (p == 0)
    {
        for (int i = 0; i < N; i++)
        {
            cout << A[i] << " ";
            if (i % n == n-1)
                cout << endl;
        }
        if (N % P != 0)
            cout << endl;
        cout << "done!" << endl;
    }

    delete [] Asub;
    delete [] A;

    // Terminate MPI
    MPI::Finalize();

    return 0;
}
