/** \file include/mpi_dispatcher/misc.hpp
** \brief Miscellaneous MPI-related definitions.
*/

#pragma once

#include <type_traits>
#include <mpi.h>

#include "../pomerol/Misc.h"

namespace pMPI {

template<typename T> MPI_Datatype mpi_datatype();

template<> MPI_Datatype mpi_datatype<Pomerol::RealType>() { return MPI_DOUBLE; }
template<> MPI_Datatype mpi_datatype<Pomerol::ComplexType>() {
  return MPI_CXX_DOUBLE_COMPLEX;
}

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
