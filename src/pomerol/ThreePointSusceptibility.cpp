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

std::ostream& operator<<(std::ostream& os, ThreePointSusceptibility::Channel channel) {
    switch(channel) {
    case ThreePointSusceptibility::PP: return os << "PP";
    case ThreePointSusceptibility::PH: return os << "PH";
    default: return os;
    }
}

auto ThreePointSusceptibility::selectChannel(QuadraticOperator const& B) -> Channel {
    auto dagger = B.getDagger();
    if(dagger == std::make_tuple(false, false))
        return PP;
    else if(dagger == std::make_tuple(true, false))
        return PH;
    else
        throw std::runtime_error(
            "ThreePointSusceptibility: Cannot select a channel based on the form of the quadratic operator B");
}

ThreePointSusceptibility::ThreePointSusceptibility(StatesClassification const& S,
                                                   Hamiltonian const& H,
                                                   CreationOperator const& CX1,
                                                   CreationOperator const& CX3,
                                                   QuadraticOperator const& Delta,
                                                   DensityMatrix const& DM)
    : Thermal(DM.beta),
      ComputableObject(),
      channel(selectChannel(Delta)),
      S(S),
      H(H),
      F1(CX1),
      F2(CX3),
      B(Delta),
      DM(DM) {
    if(channel != PP)
        throw std::runtime_error("Pairing operator for ThreePointSusceptibility must be a product of two annihilators");
}

ThreePointSusceptibility::ThreePointSusceptibility(StatesClassification const& S,
                                                   Hamiltonian const& H,
                                                   CreationOperator const& CX1,
                                                   AnnihilationOperator const& C2,
                                                   QuadraticOperator const& N,
                                                   DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), channel(selectChannel(N)), S(S), H(H), F1(CX1), F2(C2), B(N), DM(DM) {
    if(channel != PH)
        throw std::runtime_error(
            "Density operator for ThreePointSusceptibility must be a product of one creator and one annihilator");
}

void ThreePointSusceptibility::prepare() {
    if(getStatus() >= Prepared)
        return;

    // Find out non-trivial blocks of B.
    MonomialOperator::BlocksBimap const& BNontrivialBlocks = B.getBlockMapping();
    // Iterate over the outermost index.
    for(auto Biter = BNontrivialBlocks.right.begin(); Biter != BNontrivialBlocks.right.end(); Biter++) {
        BlockNumber Bright = Biter->first;
        BlockNumber Bleft = Biter->second;

        // F_1, F_2, B term
        // <Bright|F_1|F1right><F2left|F_2|Bleft><Bleft|B|Bright>
        BlockNumber F1right = F1.getRightIndex(Bright);
        BlockNumber F2left = F2.getLeftIndex(Bleft);
        if(F1right == F2left && F1right != INVALID_BLOCK_NUMBER) {
            if(DM.isRetained(Bright) && DM.isRetained(F1right) && DM.isRetained(Bleft)) {
                parts.emplace_back(F1.getPartFromLeftIndex(Bright),
                                   F2.getPartFromLeftIndex(F2left),
                                   B.getPartFromLeftIndex(Bleft),
                                   H.getPart(Bright),
                                   H.getPart(F1right),
                                   H.getPart(Bleft),
                                   DM.getPart(Bright),
                                   DM.getPart(F1right),
                                   DM.getPart(Bleft),
                                   channel == PP,
                                   false);
                parts.back().ReduceResonanceTolerance = ReduceResonanceTolerance;
            }
        }

        // F_2, F_1, B term
        // <Bright|F_2|F2right><F1left|F_1|Bleft><Bleft|B|Bright>
        BlockNumber F2right = F2.getRightIndex(Bright);
        BlockNumber F1left = F1.getLeftIndex(Bleft);
        if(F2right == F1left && F2right != INVALID_BLOCK_NUMBER) {
            if(DM.isRetained(Bright) && DM.isRetained(F2right) && DM.isRetained(Bleft)) {
                parts.emplace_back(F2.getPartFromLeftIndex(Bright),
                                   F1.getPartFromLeftIndex(F1left),
                                   B.getPartFromLeftIndex(Bleft),
                                   H.getPart(Bright),
                                   H.getPart(F2right),
                                   H.getPart(Bleft),
                                   DM.getPart(Bright),
                                   DM.getPart(F2right),
                                   DM.getPart(Bleft),
                                   channel == PP,
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
    switch(channel) {
    case PP: {
        switch(Position) {
        case 0: return F1.getIndex();
        case 1: return B.getIndex1();
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
        case 2: return static_cast<CreationOperator const&>(F2).getIndex();
        case 3: return B.getIndex2();
        }
    }
    case PH: {
        switch(Position) {
        case 0: return F1.getIndex();
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
        case 1: return static_cast<AnnihilationOperator const&>(F2).getIndex();
        case 2: return B.getIndex1();
        case 3: return B.getIndex2();
        }
    }
    default: assert(0);
    }
    throw std::runtime_error("ThreePointSusceptibility: Could not get operator index");
}

} // namespace Pomerol
