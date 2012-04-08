#include "config.h"
#include <iostream>

#if defined(USING_MPI)
#include "mpi.h"

int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);

    try {
        int rank = MPI::COMM_WORLD.Get_rank();
        int size = MPI::COMM_WORLD.Get_size();
        std::cout << "Hello world from process " << rank << " of " << size << std::endl;
    } catch (MPI::Exception e) {
        std::cout << "MPI ERROR: " << e.Get_error_code()
                  << " - " << e.Get_error_string()
                  << std::endl;
    }

    MPI::Finalize();
    return 0;
}

#else

int main(void)
{
    std::cerr << "MPI not available" << std::endl;
    return 1;
}

#endif // USING_MPI