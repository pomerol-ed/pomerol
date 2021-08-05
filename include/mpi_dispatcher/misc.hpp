/** \file include/mpi_dispatcher/misc.hpp
** \brief Miscellaneous MPI-related definitions.
*/
#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MISC_H
#define POMEROL_INCLUDE_MPI_DISPATCHER_MISC_H

#include <mpi.h>

namespace pMPI {

// Size of communicator 'Comm'
inline int size(MPI_Comm const& Comm) {
    int s;
    MPI_Comm_size(Comm, &s);
    return s;
}

// Rank of this process in communicator 'Comm'
inline int rank(MPI_Comm const& Comm) {
    int r;
    MPI_Comm_rank(Comm, &r);
    return r;
}

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MISC_H
