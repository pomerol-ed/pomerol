//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/TwoParticleGF.cpp
/// \brief Fermionic two-particle Matsubara Green's function (implementation).
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/TwoParticleGF.hpp"

#include "mpi_dispatcher/mpi_skel.hpp"

#include <array>
#include <cassert>
#include <map>
#include <stdexcept>

namespace Pomerol {

TwoParticleGF::TwoParticleGF(StatesClassification const& S,
                             Hamiltonian const& H,
                             AnnihilationOperator const& C1,
                             AnnihilationOperator const& C2,
                             CreationOperator const& CX3,
                             CreationOperator const& CX4,
                             DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM) {}

BlockNumber
TwoParticleGF::getLeftIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber RightIndex) const {
    switch(permutations3[PermutationNumber].perm[OperatorPosition]) {
    case 0: return C1.getLeftIndex(RightIndex);
    case 1: return C2.getLeftIndex(RightIndex);
    case 2: return CX3.getLeftIndex(RightIndex);
    default: return INVALID_BLOCK_NUMBER;
    }
}

BlockNumber
TwoParticleGF::getRightIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber LeftIndex) const {
    switch(permutations3[PermutationNumber].perm[OperatorPosition]) {
    case 0: return C1.getRightIndex(LeftIndex);
    case 1: return C2.getRightIndex(LeftIndex);
    case 2: return CX3.getRightIndex(LeftIndex);
    default: return INVALID_BLOCK_NUMBER;
    }
}

MonomialOperatorPart const& TwoParticleGF::OperatorPartAtPosition(std::size_t PermutationNumber,
                                                                  std::size_t OperatorPosition,
                                                                  BlockNumber LeftIndex) const {
    switch(permutations3[PermutationNumber].perm[OperatorPosition]) {
    case 0: return C1.getPartFromLeftIndex(LeftIndex);
    case 1: return C2.getPartFromLeftIndex(LeftIndex);
    case 2: return CX3.getPartFromLeftIndex(LeftIndex);
    default: throw std::runtime_error("TwoParticleGF: Could not find operator part");
    }
}

void TwoParticleGF::prepare() {
    if(getStatus() >= Prepared)
        return;

    // Find out non-trivial blocks of CX4.
    MonomialOperator::BlocksBimap const& CX4NontrivialBlocks = CX4.getBlockMapping();
    for(auto outer_iter = CX4NontrivialBlocks.right.begin(); outer_iter != CX4NontrivialBlocks.right.end();
        outer_iter++) {                      // Iterate over the outermost index.
        for(std::size_t p = 0; p < 6; ++p) { // Choose a permutation
            std::array<BlockNumber, 4> LeftIndices{};
            LeftIndices[0] = outer_iter->first;
            LeftIndices[3] = outer_iter->second;
            LeftIndices[2] = getLeftIndex(p, 2, LeftIndices[3]);
            LeftIndices[1] = getRightIndex(p, 0, LeftIndices[0]);
            // < LeftIndices[0] | O_1 | LeftIndices[1] >
            // < LeftIndices[1] | O_2 | getRightIndex(p,1,LeftIndices[1]) >
            // < LeftIndices[2]| O_3 | LeftIndices[3] >
            // < LeftIndices[3] | CX4 | LeftIndices[0] >
            // Select a relevant 'world stripe' (sequence of blocks).
            if(getRightIndex(p, 1, LeftIndices[1]) == LeftIndices[2] && (LeftIndices[1] != INVALID_BLOCK_NUMBER) &&
               (LeftIndices[2] != INVALID_BLOCK_NUMBER)) {
                { // check if retained blocks are included. If not, do not push.
                    bool include_block_retained = false;
                    for(int k = 0; k < 4; ++k)
                        if(DM.isRetained(LeftIndices[k]))
                            include_block_retained = true;
                    if(!include_block_retained)
                        continue;
                }
                parts.emplace_back(OperatorPartAtPosition(p, 0, LeftIndices[0]),
                                   OperatorPartAtPosition(p, 1, LeftIndices[1]),
                                   OperatorPartAtPosition(p, 2, LeftIndices[2]),
                                   (MonomialOperatorPart&)CX4.getPartFromLeftIndex(LeftIndices[3]),
                                   H.getPart(LeftIndices[0]),
                                   H.getPart(LeftIndices[1]),
                                   H.getPart(LeftIndices[2]),
                                   H.getPart(LeftIndices[3]),
                                   DM.getPart(LeftIndices[0]),
                                   DM.getPart(LeftIndices[1]),
                                   DM.getPart(LeftIndices[2]),
                                   DM.getPart(LeftIndices[3]),
                                   permutations3[p]);

                parts.back().ReduceResonanceTolerance = ReduceResonanceTolerance;
                parts.back().CoefficientTolerance = CoefficientTolerance;
                parts.back().MultiTermCoefficientTolerance = MultiTermCoefficientTolerance;
            }
        }
    }
    if(!parts.empty()) {
        Vanishing = false;
        INFO("TwoParticleGF(" << getIndex(0) << getIndex(1) << getIndex(2) << getIndex(3) << "): " << parts.size()
                              << " parts will be calculated");
    }
    setStatus(Prepared);
}

// An mpi adapter to 1) compute 2pgf terms; 2) convert them to a Matsubara Container; 3) purge terms
struct ComputeAndClearWrap {
    ComputeAndClearWrap(FreqVec const& freqs,
                        std::vector<ComplexType>& data,
                        TwoParticleGFPart& p,
                        bool clear,
                        bool fill,
                        int complexity = 1)
        : complexity(complexity), freqs_(freqs), data_(data), p(p), clear_(clear), fill_(fill) {}

    void run() {
        p.compute();
        if(fill_) {
            std::size_t wsize = freqs_.size();
#ifdef POMEROL_USE_OPENMP
#pragma omp parallel for
#endif
            for(int w = 0; w < wsize; ++w) {
                data_[w] += p(std::get<0>(freqs_[w]), std::get<1>(freqs_[w]), std::get<2>(freqs_[w]));
            }
#ifdef POMEROL_USE_OPENMP
#pragma omp barrier
#endif
        }
        if(clear_)
            p.clear();
    }

    // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
    int const complexity;

private:
    FreqVec const& freqs_;
    std::vector<ComplexType>& data_;
    TwoParticleGFPart& p;
    bool clear_;
    bool fill_;
};

std::vector<ComplexType> TwoParticleGF::compute(bool clear, FreqVec const& freqs, MPI_Comm const& comm) {
    if(getStatus() < Prepared)
        throw StatusMismatch("TwoParticleGF is not prepared yet.");

    std::vector<ComplexType> m_data;
    if(getStatus() >= Computed)
        return m_data;

    if(!Vanishing) {
        // Create a "skeleton" class with pointers to part that can call a compute method
        pMPI::mpi_skel<ComputeAndClearWrap> skel;
        bool fill_container = !freqs.empty();
        skel.parts.reserve(parts.size());
        m_data.resize(freqs.size(), 0.0);
        for(auto& part : parts) {
            skel.parts.emplace_back(freqs, m_data, part, clear, fill_container, 1);
        }
        std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true); // actual running - very costly

        // Start distributing data
        MPI_Barrier(comm);

        MPI_Allreduce(MPI_IN_PLACE,
                      m_data.data(),
                      static_cast<int>(m_data.size()),
                      MPI_CXX_DOUBLE_COMPLEX,
                      MPI_SUM,
                      comm);

        // Optionally distribute terms to other processes
        if(!clear) {
            for(int p = 0; p < static_cast<int>(parts.size()); ++p) {
                parts[p].NonResonantTerms.broadcast(comm, job_map[p]);
                parts[p].ResonantTerms.broadcast(comm, job_map[p]);
                parts[p].setStatus(TwoParticleGFPart::Computed);
            }
            MPI_Barrier(comm);
        }
    }

    setStatus(Computed);

    return m_data;
}

ParticleIndex TwoParticleGF::getIndex(std::size_t Position) const {
    switch(Position) {
    case 0: return C1.getIndex();
    case 1: return C2.getIndex();
    case 2: return CX3.getIndex();
    case 3: return CX4.getIndex();
    default: assert(0);
    }
    throw std::runtime_error("TwoParticleGF: Could not get operator index");
}

} // namespace Pomerol
