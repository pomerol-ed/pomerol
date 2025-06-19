//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/GreensFunction.cpp
/// \brief Fermionic single-particle Matsubara Green's function (implementation).
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/GreensFunction.hpp"

#include <cassert>
#include <stdexcept>

namespace Pomerol {

GreensFunction::GreensFunction(StatesClassification const& S,
                               Hamiltonian const& H,
                               FieldOperator const& F1,
                               FieldOperator const& F2,
                               DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), S(S), H(H), F1(F1), F2(F2), DM(DM) {}

GreensFunction::GreensFunction(GreensFunction const& GF)
    : Thermal(GF.beta),
      ComputableObject(GF),
      S(GF.S),
      H(GF.H),
      F1(GF.F1),
      F2(GF.F2),
      DM(GF.DM),
      Vanishing(GF.Vanishing),
      parts(GF.parts) {}

void GreensFunction::prepare() {
    if(getStatus() >= Prepared)
        return;

    // Find out non-trivial blocks of F1 and F2.
    MonomialOperator::BlocksBimap const& F1NontrivialBlocks = F1.getBlockMapping();
    MonomialOperator::BlocksBimap const& F2NontrivialBlocks = F2.getBlockMapping();

    auto F1iter = F1NontrivialBlocks.left.begin();
    auto F2iter = F2NontrivialBlocks.right.begin();

    while(F1iter != F1NontrivialBlocks.left.end() && F2iter != F2NontrivialBlocks.right.end()) {
        // <F1left|F1|F1right><F2left|F2|F2right>
        BlockNumber F1left = F1iter->first;
        BlockNumber F1right = F1iter->second;
        BlockNumber F2left = F2iter->second;
        BlockNumber F2right = F2iter->first;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(F1left == F2right && F1right == F2left) {
            // check if retained blocks are included. If not, do not push.
            if(DM.isRetained(F1left) || DM.isRetained(F1right)) {
                parts.emplace_back(F1.getPartFromLeftIndex(F1left),
                                   F2.getPartFromRightIndex(F2right),
                                   H.getPart(F1right),
                                   H.getPart(F1left),
                                   DM.getPart(F1right),
                                   DM.getPart(F1left),
                                   PoleResolution,
                                   CoefficientTolerance);
            }
        }

        unsigned long F1leftInt = F1left;
        unsigned long F2rightInt = F2right;

        if(F1leftInt <= F2rightInt)
            F1iter++;
        if(F1leftInt >= F2rightInt)
            F2iter++;
    }
    if(!parts.empty())
        Vanishing = false;

    setStatus(Prepared);
}

void GreensFunction::compute() {
    if(getStatus() >= Computed)
        return;
    if(getStatus() < Prepared)
        prepare();

    if(getStatus() < Computed) {
        for(auto& p : parts)
            p.compute();
    }

    setStatus(Computed);
}

ParticleIndex GreensFunction::getIndex(std::size_t Position) const {
    switch(Position) {
    case 0: return F1.getIndex();
    case 1: return F2.getIndex();
    default: throw std::runtime_error("GreensFunction: Wrong operator");
    }
}

} // namespace Pomerol
