//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/ThreePointSusceptibility.cpp
/// \brief 3-point susceptibility in the Matsubara representation (implementation).
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include "pomerol/ThreePointSusceptibility.hpp"

#include "mpi_dispatcher/mpi_skel.hpp"

#include <array>
#include <cassert>
#include <ostream>
#include <stdexcept>
#include <tuple>

namespace Pomerol {

ThreePointSusceptibility::ThreePointSusceptibility(Channel channel,
                                                   StatesClassification const& S,
                                                   Hamiltonian const& H,
                                                   CreationOperator const& CX1,
                                                   AnnihilationOperator const& C2,
                                                   CreationOperator const& CX3,
                                                   AnnihilationOperator const& C4,
                                                   DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), channel(channel), S(S), H(H), CX1(CX1), C2(C2), CX3(CX3), C4(C4), DM(DM) {}

CreationOperator const& ThreePointSusceptibility::getF1() const {
    return CX1;
}
MonomialOperator const& ThreePointSusceptibility::getF2() const {
    switch(channel) {
    case PP: return CX3;
    case PH: return C2;
    case xPH: return C4;
    }
}
MonomialOperator const& ThreePointSusceptibility::getB1() const {
    switch(channel) {
    case PP: return C2;
    case PH: return CX3;
    case xPH: return CX3;
    }
}
MonomialOperator const& ThreePointSusceptibility::getB2() const {
    switch(channel) {
    case PP: return C4;
    case PH: return C4;
    case xPH: return C2;
    }
}

void ThreePointSusceptibility::prepare() {
    if(getStatus() >= Prepared)
        return;

    CreationOperator const& F1 = getF1();
    MonomialOperator const& F2 = getF2();
    MonomialOperator const& B1 = getB1();
    MonomialOperator const& B2 = getB2();

    // Find out non-trivial blocks of B2.
    MonomialOperator::BlocksBimap const& B2NontrivialBlocks = B2.getBlockMapping();

    // Iterate over the outermost index.
    for(auto B2iter = B2NontrivialBlocks.right.begin(); B2iter != B2NontrivialBlocks.right.end(); B2iter++) {
        BlockNumber B2right = B2iter->first;
        BlockNumber B2left = B2iter->second;
        BlockNumber B1left = B1.getLeftIndex(B2left);
        if(B1left == INVALID_BLOCK_NUMBER)
            continue;

        // F_1, F_2, B term
        // <B2right|F_1|F1right><F2left|F_2|B1left><B1left|B_1|B2left><B2left|B_2|B2right>
        BlockNumber F1right = F1.getRightIndex(B2right);
        BlockNumber F2left = F2.getLeftIndex(B1left);
        if(F1right == F2left && F1right != INVALID_BLOCK_NUMBER) {
            if(DM.isRetained(B2right) && DM.isRetained(F1right) && DM.isRetained(B1left)) {
                parts.emplace_back(F1.getPartFromLeftIndex(B2right),
                                   F2.getPartFromLeftIndex(F2left),
                                   B1.getPartFromLeftIndex(B1left),
                                   B2.getPartFromLeftIndex(B2left),
                                   H.getPart(B2right),
                                   H.getPart(F1right),
                                   H.getPart(B1left),
                                   DM.getPart(B2right),
                                   DM.getPart(F1right),
                                   DM.getPart(B1left),
                                   channel,
                                   false);
                parts.back().ReduceResonanceTolerance = ReduceResonanceTolerance;
            }
        }

        // F_2, F_1, B term
        // <B2right|F_2|F2right><F1left|F_1|B1left><B1left|B_1|B2left><B2left|B_2|B2right>
        BlockNumber F2right = F2.getRightIndex(B2right);
        BlockNumber F1left = F1.getLeftIndex(B1left);

        if(F2right == F1left && F2right != INVALID_BLOCK_NUMBER) {
            if(DM.isRetained(B2right) && DM.isRetained(F2right) && DM.isRetained(B1left)) {
                parts.emplace_back(F2.getPartFromLeftIndex(B2right),
                                   F1.getPartFromLeftIndex(F1left),
                                   B1.getPartFromLeftIndex(B1left),
                                   B2.getPartFromLeftIndex(B2left),
                                   H.getPart(B2right),
                                   H.getPart(F2right),
                                   H.getPart(B1left),
                                   DM.getPart(B2right),
                                   DM.getPart(F2right),
                                   DM.getPart(B1left),
                                   channel,
                                   true);
                parts.back().ReduceResonanceTolerance = ReduceResonanceTolerance;
                parts.back().CoefficientTolerance = CoefficientTolerance;
            }
        }
    }

    if(!parts.empty()) {
        Vanishing = false;
        INFO("ThreePointSusceptibility(" << channel << "," << getIndex(0) << "," << getIndex(1) << "," << getIndex(2)
                                         << "," << getIndex(3) << "): " << parts.size() << " parts will be calculated");
    }
    setStatus(Prepared);
}

// An MPI adapter to
// 1) compute ThreePointSusceptibility terms;
// 2) convert them to a Matsubara Container;
// 3) purge terms
struct ComputeAndClearWrap3PSusc {
    ComputeAndClearWrap3PSusc(FreqVec2 const& freqs,
                              std::vector<ComplexType>& data,
                              ThreePointSusceptibilityPart& p,
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
                data_[w] += p(std::get<0>(freqs_[w]), std::get<1>(freqs_[w]));
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
    FreqVec2 const& freqs_;
    std::vector<ComplexType>& data_;
    ThreePointSusceptibilityPart& p;
    bool clear_;
    bool fill_;
};

std::vector<ComplexType> ThreePointSusceptibility::compute(bool clear, FreqVec2 const& freqs, MPI_Comm const& comm) {
    if(getStatus() < Prepared)
        throw StatusMismatch("ThreePointSusceptibility is not prepared yet.");

    std::vector<ComplexType> m_data;
    if(getStatus() >= Computed)
        return m_data;

    if(!Vanishing) {
        // Create a "skeleton" class with pointers to part that can call a compute method
        pMPI::mpi_skel<ComputeAndClearWrap3PSusc> skel;
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
                      POMEROL_MPI_DOUBLE_COMPLEX,
                      MPI_SUM,
                      comm);

        // Optionally distribute terms to other processes
        if(!clear) {
            for(int p = 0; p < static_cast<int>(parts.size()); ++p) {
                parts[p].NonResonantFFTerms.broadcast(comm, job_map[p]);
                parts[p].NonResonantFBTerms.broadcast(comm, job_map[p]);
                parts[p].ResonantTerms.broadcast(comm, job_map[p]);
                parts[p].setStatus(ThreePointSusceptibilityPart::Computed);
            }
            MPI_Barrier(comm);
        }
    }

    setStatus(Computed);

    return m_data;
}

ParticleIndex ThreePointSusceptibility::getIndex(std::size_t Position) const {
    switch(Position) {
    case 0: return CX1.getIndex();
    case 1: return C2.getIndex();
    case 2: return CX3.getIndex();
    case 3: return C4.getIndex();
    default: assert(0);
    }
    throw std::runtime_error("ThreePointSusceptibility: Could not get operator index");
}

} // namespace Pomerol
