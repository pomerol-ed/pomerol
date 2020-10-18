/** \file include/mpi_dispatcher/misc.hpp
** \brief Miscellaneous MPI-related definitions.
*/

#pragma once

#include <mpi.h>

#include "../pomerol/Misc.h"

#ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
#define MPI_MELEM_DATATYPE MPI_CXX_DOUBLE_COMPLEX
#else
#define MPI_MELEM_DATATYPE MPI_DOUBLE
#endif

namespace pMPI {

// Size of communicator 'Comm'
inline int size(const MPI_Comm &Comm) {
    int s;
    MPI_Comm_size(Comm, &s);
    return s;
}

// Rank of this process in communicator 'Comm'
inline int rank(const MPI_Comm &Comm) {
    int r;
    MPI_Comm_rank(Comm, &r);
    return r;
}

} // end of namespace pMPI
