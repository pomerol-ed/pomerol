//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/FieldOperatorContainer.cpp
/// \brief A container for creation and annihilation operators (implementation).
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/FieldOperatorContainer.hpp"

#include <stdexcept>

namespace Pomerol {

void FieldOperatorContainer::computeAll() {
    for(auto& cdag_p : mapCreationOperators) {
        auto& cdag = cdag_p.second;
        cdag.compute();

        auto& c = mapAnnihilationOperators.find(cdag_p.first)->second;

        auto const& cdag_block_map = cdag_p.second.getBlockMapping();
        for(auto cdag_map_it = cdag_block_map.right.begin(); cdag_map_it != cdag_block_map.right.end(); ++cdag_map_it) {
            auto& cPart = c.getPartFromRightIndex(cdag_map_it->second);
            auto& cdagPart = cdag.getPartFromRightIndex(cdag_map_it->first);
            cPart.setFromAdjoint(cdagPart);
        }
        c.setStatus(ComputableObject::Computed);
    }
}

CreationOperator const& FieldOperatorContainer::getCreationOperator(ParticleIndex in) const {
    auto it = mapCreationOperators.find(in);
    if(it == mapCreationOperators.end())
        throw std::runtime_error("Creation operator not found.");
    else
        return it->second;
}

AnnihilationOperator const& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in) const {
    auto it = mapAnnihilationOperators.find(in);
    if(it == mapAnnihilationOperators.end())
        throw std::runtime_error("Annihilation operator not found.");
    else
        return it->second;
}

} // namespace Pomerol
