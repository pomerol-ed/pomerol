//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/Hamiltonian.cpp
/// \brief Storage and diagonalization of a Hamiltonian matrix (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include "pomerol/Hamiltonian.hpp"

#include "mpi_dispatcher/mpi_skel.hpp"

#include <cstddef>
#include <map>
#include <stdexcept>

namespace Pomerol {

// cppcheck-suppress unusedPrivateFunction
template <bool C> void Hamiltonian::prepareImpl(LOperatorTypeRC<C> const& HOp, MPI_Comm const& comm) {
    BlockNumber NumberOfBlocks = S.getNumberOfBlocks();
    int comm_rank = pMPI::rank(comm);
    if(!comm_rank)
        INFO_NONEWLINE("Preparing Hamiltonian parts...");

    parts.reserve(NumberOfBlocks);
    for(BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; ++CurrentBlock) {
        parts.emplace_back(HOp, S, CurrentBlock);
    }

    pMPI::mpi_skel<pMPI::PrepareWrap<HamiltonianPart>> skel;
    skel.parts.reserve(parts.size());
    for(auto& part : parts) {
        skel.parts.emplace_back(part);
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, false);
    MPI_Barrier(comm);

    MPI_Datatype H_dt = C ? POMEROL_MPI_DOUBLE_COMPLEX : MPI_DOUBLE;

    for(int p = 0; p < static_cast<int>(parts.size()); ++p) {
        auto& part = parts[p];
        if(comm_rank == job_map[p]) {
            if(part.getStatus() != HamiltonianPart::Prepared) {
                ERROR("Worker" << comm_rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            auto& H = part.getMatrix<C>();
            MPI_Bcast(H.data(), H.size(), H_dt, comm_rank, comm);
        } else {
            part.initHMatrix<C>();
            auto& H = part.getMatrix<C>();
            MPI_Bcast(H.data(), H.rows() * H.cols(), H_dt, job_map[p], comm);
            part.setStatus(HamiltonianPart::Prepared);
        }
    }
}

template void Hamiltonian::prepareImpl<true>(LOperatorTypeRC<true> const&, MPI_Comm const&);
template void Hamiltonian::prepareImpl<false>(LOperatorTypeRC<false> const&, MPI_Comm const&);

template <bool C> void Hamiltonian::computeImpl(MPI_Comm const& comm) {
    // Create a "skeleton" class with pointers to part that can call a compute method
    pMPI::mpi_skel<pMPI::ComputeWrap<HamiltonianPart>> skel;
    skel.parts.reserve(parts.size());
    for(auto& part : parts) {
        skel.parts.emplace_back(part, static_cast<int>(part.getSize()));
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true);
    int comm_rank = pMPI::rank(comm);

    // Start distributing data
    MPI_Barrier(comm);
    MPI_Datatype H_dt = C ? POMEROL_MPI_DOUBLE_COMPLEX : MPI_DOUBLE;
    for(int p = 0; p < static_cast<int>(parts.size()); ++p) {
        auto& part = parts[p];
        auto& H = part.getMatrix<C>();
        if(comm_rank == job_map[p]) {
            if(part.getStatus() != HamiltonianPart::Computed) {
                ERROR("Worker" << comm_rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            MPI_Bcast(H.data(), H.size(), H_dt, comm_rank, comm);
            MPI_Bcast(part.Eigenvalues.data(), static_cast<int>(part.Eigenvalues.size()), MPI_DOUBLE, comm_rank, comm);
        } else {
            part.Eigenvalues.resize(H.rows());
            MPI_Bcast(H.data(), H.size(), H_dt, job_map[p], comm);
            MPI_Bcast(part.Eigenvalues.data(), static_cast<int>(part.Eigenvalues.size()), MPI_DOUBLE, job_map[p], comm);
            part.setStatus(HamiltonianPart::Computed);
        }
    }
}

void Hamiltonian::compute(MPI_Comm const& comm) {
    if(getStatus() >= Computed)
        return;

    if(Complex)
        computeImpl<true>(comm);
    else
        computeImpl<false>(comm);

    computeGroundEnergy();

    setStatus(Computed);
}

void Hamiltonian::reduce(RealType Cutoff) {
    INFO("Performing EV cutoff at " << Cutoff << " level");
    for(auto& part : parts)
        part.reduce(GroundEnergy + Cutoff);
}

InnerQuantumState Hamiltonian::getBlockSize(BlockNumber Block) const {
    return parts[Block].getSize();
}

void Hamiltonian::computeGroundEnergy() {
    RealVectorType lev(S.getNumberOfBlocks());
    for(BlockNumber b = 0; b < parts.size(); ++b)
        lev(b) = parts[b].getMinimumEigenvalue();
    GroundEnergy = lev.minCoeff();
}

RealType Hamiltonian::getEigenValue(QuantumState state) const {
    return parts[S.getBlockNumber(state)].getEigenValue(S.getInnerState(state));
}

RealVectorType const& Hamiltonian::getEigenValues(BlockNumber Block) const {
    return parts[Block].getEigenValues();
}

RealVectorType Hamiltonian::getEigenValues() const {
    RealVectorType out(S.getNumberOfStates());
    long copied_size = 0;
    for(auto const& part : parts) {
        auto const& ev = part.getEigenValues();
        out.segment(copied_size, ev.size()) = ev;
        copied_size += ev.size();
    }
    return out;
}

} // namespace Pomerol
