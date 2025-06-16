//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/ThreePointSusceptibilityContainer.cpp
/// \brief Storage for multiple 3-point susceptibilities in the Matsubara representation (implementation).
/// \author Igor Krivenko

#include "pomerol/ThreePointSusceptibilityContainer.hpp"

namespace Pomerol {

void ThreePointSusceptibilityContainer::prepareAll(std::set<IndexCombination4> const& InitialIndices) {
    fill(InitialIndices);
    for(auto& el : ElementsMap) {
        auto& chi3 = static_cast<ThreePointSusceptibility&>(el.second);
        chi3.PoleResolution = PoleResolution;
        chi3.CoefficientTolerance = CoefficientTolerance;
        chi3.prepare();
    }
}

std::map<IndexCombination4, std::vector<ComplexType>>
ThreePointSusceptibilityContainer::computeAll(bool clearTerms, FreqVec2 const& freqs, MPI_Comm const& comm) {
    std::map<IndexCombination4, std::vector<ComplexType>> out;
    for(auto& el : ElementsMap) {
        INFO("Computing 3PSusceptibility for " << el.first);
        auto& chi3 = static_cast<ThreePointSusceptibility&>(el.second);
        out.emplace(el.first, chi3.compute(clearTerms, freqs, comm));
    }
    return out;
}

std::shared_ptr<ThreePointSusceptibility>
ThreePointSusceptibilityContainer::createElement(IndexCombination4 const& Indices) const {
    CreationOperator const& CX1 = Operators.getCreationOperator(Indices.Index1);
    AnnihilationOperator const& C2 = Operators.getAnnihilationOperator(Indices.Index2);
    CreationOperator const& CX3 = Operators.getCreationOperator(Indices.Index3);
    AnnihilationOperator const& C4 = Operators.getAnnihilationOperator(Indices.Index4);

    return std::make_shared<ThreePointSusceptibility>(channel, S, H, CX1, C2, CX3, C4, DM);
}

} // namespace Pomerol
