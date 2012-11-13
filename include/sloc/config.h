#ifndef SLOC_CONFIG_H
#define SLOC_CONFIG_H

#include <deal.II/base/config.h>

// check for MPI support
#if defined(DEAL_II_COMPILER_SUPPORTS_MPI)
#define HAVE_MPI_H  1
#define USING_MPI   1
#endif

#endif
