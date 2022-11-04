//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/mpi_dispatcher/misc.hpp
/// \brief Miscellaneous MPI-related functions.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MISC_HPP
#define POMEROL_INCLUDE_MPI_DISPATCHER_MISC_HPP

#include <mpi.h>

/// A bunch of tools used for MPI-parallelization of computations.
namespace pMPI {

/// \defgroup MPI MPI tools
///@{

/// Return size of an MPI communicator.
/// \param[in] Comm MPI communicator.
inline int size(MPI_Comm const& Comm) {
    int s;
    MPI_Comm_size(Comm, &s);
    return s;
}

/// Return the rank of the calling process in a given MPI communicator.
/// \param[in] Comm MPI communicator.
inline int rank(MPI_Comm const& Comm) {
    int r;
    MPI_Comm_rank(Comm, &r);
    return r;
}

///@}

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MISC_HPP
