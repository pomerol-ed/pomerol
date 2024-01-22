//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/StatesClassification.cpp
/// \brief Classification of many-body basis states (Fock states) into subspaces (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include "pomerol/StatesClassification.hpp"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>

namespace Pomerol {

//
// StatesClassification
//

void StatesClassification::initSingleBlock(QuantumState Dim) {
    StateBlockIndex.resize(Dim, 0);
    StatesContainer.emplace_back(Dim, 0);
    std::iota(StatesContainer.back().begin(), StatesContainer.back().end(), 0);
}

void StatesClassification::initMultipleBlocks(libcommute::space_partition const& partition) {
    StateBlockIndex.resize(partition.dim());
    StatesContainer.resize(partition.n_subspaces(), std::vector<QuantumState>());
    foreach(partition, [this](QuantumState State, BlockNumber Block) {
        StateBlockIndex[State] = Block;
        StatesContainer[Block].push_back(State);
    })
        ;
}

void StatesClassification::checkComputed() const {
    if(getStatus() < Computed) {
        throw StatusMismatch("StatesClassification is not computed yet.");
    }
}

Pomerol::InnerQuantumState StatesClassification::getBlockSize(BlockNumber in) const {
    checkComputed();
    return static_cast<InnerQuantumState>(getFockStates(in).size());
}

std::vector<QuantumState> const& StatesClassification::getFockStates(BlockNumber in) const {
    checkComputed();
    return StatesContainer[in];
}

QuantumState StatesClassification::getFockState(BlockNumber in, InnerQuantumState m) const {
    checkComputed();
    if(int(in) < StatesContainer.size())
        if(m < StatesContainer[in].size())
            return StatesContainer[in][m];
    throw std::runtime_error("Wrong inner state " + std::to_string(m));
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const {
    checkComputed();
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    return StateBlockIndex[in];
}

InnerQuantumState StatesClassification::getInnerState(QuantumState in) const {
    checkComputed();
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    BlockNumber n = StateBlockIndex[in];
    auto const& Block = StatesContainer[n];
    return std::distance(Block.begin(), std::find(Block.begin(), Block.end(), in));
}

} // namespace Pomerol
