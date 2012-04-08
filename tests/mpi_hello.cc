#include "config.h"
#include <stdio.h>
#if defined(USING_MPI)

#include "mpi.h"
int main(int argc, char *argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
	MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
	printf("Hello world from process %d of %d\n", rank, size);
	MPI_Finalize();
	return 0;
}

#else

int main(void)
{
    fprintf(stderr, "MPI not available\n");
    return 1;
}

#endif // USING_MPI
