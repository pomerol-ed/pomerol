//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/GFContainer.cpp
/// \brief Storage for multiple fermionic single-particle Matsubara Green's functions (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include "pomerol/GFContainer.hpp"

namespace Pomerol {

void GFContainer::prepareAll(std::set<IndexCombination2> const& Indices) {
    fill(Indices);
    for(auto& el : ElementsMap) {
        el.second->MatrixElementTolerance = MatrixElementTolerance;
        el.second->prepare();
    }
}

void GFContainer::computeAll() {
    for(auto& el : ElementsMap)
        el.second->compute();
}

std::shared_ptr<GreensFunction> GFContainer::createElement(IndexCombination2 const& Indices) const {
    return std::make_shared<GreensFunction>(S,
                                            H,
                                            Operators.getAnnihilationOperator(Indices.Index1),
                                            Operators.getCreationOperator(Indices.Index2),
                                            DM);
}

} // namespace Pomerol
