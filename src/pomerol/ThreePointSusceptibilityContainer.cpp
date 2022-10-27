//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/ThreePointSusceptibilityContainer.cpp
/// \brief Storage for multiple 3-point susceptibilities in the Matsubara representation (implementation).
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include "pomerol/ThreePointSusceptibilityContainer.hpp"

namespace Pomerol {

void ThreePointSusceptibilityContainer::prepareAll(std::set<IndexCombination2> const& InitialIndices) {
    fill(InitialIndices);
    for(auto& el : ElementsMap) {
        el.second->ReduceResonanceTolerance = ReduceResonanceTolerance;
        el.second->CoefficientTolerance = CoefficientTolerance;
        el.second->prepare();
    }
}

std::map<IndexCombination2, std::vector<ComplexType>>
ThreePointSusceptibilityContainer::computeAll(bool clear, FreqVec2 const& freqs, MPI_Comm const& comm) {
    std::map<IndexCombination2, std::vector<ComplexType>> out;
    for(auto& el : ElementsMap) {
        INFO("Computing 3PSusceptibility for " << el.first);
        out.emplace(el.first, el.second->compute(clear, freqs, comm));
    }
    return out;
}

std::shared_ptr<ThreePointSusceptibility>
ThreePointSusceptibilityContainer::createElement(IndexCombination2 const& Indices) const {
    CreationOperator const& CX1 = Operators.getCreationOperator(Indices.Index1);
    switch(Channel) {
    case ThreePointSusceptibility::PP:
        return std::make_shared<ThreePointSusceptibility>(S,
                                                          H,
                                                          CX1,
                                                          Operators.getCreationOperator(Indices.Index2),
                                                          B,
                                                          DM);
    case ThreePointSusceptibility::PH:
    case ThreePointSusceptibility::CrossedPH:
        return std::make_shared<ThreePointSusceptibility>(S,
                                                          H,
                                                          CX1,
                                                          Operators.getAnnihilationOperator(Indices.Index2),
                                                          B,
                                                          DM);
    }
}

} // namespace Pomerol
