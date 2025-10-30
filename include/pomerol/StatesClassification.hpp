//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/StatesClassification.hpp
/// \brief Classification of many-body basis states (Fock states) into subspaces.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_STATESCLASSIFICATION_HPP
#define POMEROL_INCLUDE_STATESCLASSIFICATION_HPP

#include "HilbertSpace.hpp"
#include "Misc.hpp"

#include <libcommute/loperator/space_partition.hpp>

#include <vector>

namespace Pomerol {

/// \addtogroup ED
///@{

/// Index of a subspace (block) within a full many-body Hilbert space.
using BlockNumber = int;
/// A special value that stands for a non-existent subspace (block).
constexpr BlockNumber INVALID_BLOCK_NUMBER = -1;

/// Index of a state within a block.
using InnerQuantumState = libcommute::sv_index_type;

/// \brief Classification of many-body basis states into bases of invariant subspaces.
///
/// This class stores lists of Fock states belonging to each invariant subspace (block)
/// of a Hilbert space.
class StatesClassification : public ComputableObject {

    /// Lists of Fock states spanning the invariant subspaces, one inner vector per subspace.
    std::vector<std::vector<QuantumState>> StatesContainer;
    /// Each element of this vector is the block number the corresponding Fock state belongs to.
    std::vector<BlockNumber> StateBlockIndex;

public:
    /// Construct without filling any Fock state lists.
    StatesClassification() = default;

    /// Populate the Fock state lists from a \ref HilbertSpace object. If the \ref HilbertSpace is
    /// not in the \ref ComputableObject::Computed state, then existence of just one invariant
    /// subspace coinciding with the full Hilbert space will be assumed.
    /// \tparam IndexTypes Types of indices carried by operators acting in the Hilbert space \p HS.
    /// \param[in] HS The Hilbert space.
    template <typename... IndexTypes> void compute(HilbertSpace<IndexTypes...> const& HS) {
        if(getStatus() == Computed)
            return;
        auto const& FullHilbertSpace = HS.getFullHilbertSpace();
        auto Dim = FullHilbertSpace.dim();
        if(HS.getStatus() == Computed) { // Multiple blocks revealed by HS
            initMultipleBlocks(HS.getSpacePartition());
        } else { // Just one block
            initSingleBlock(Dim);
        }
        setStatus(Computed);
    }

    /// Get the total number of Fock states.
    QuantumState getNumberOfStates() const { return StateBlockIndex.size(); }

    /// Get the number of the invariant subspaces.
    BlockNumber getNumberOfBlocks() const { return StatesContainer.size(); }

    /// Get the number of Fock states spanning a given invariant subspace.
    /// \param[in] in Index of the invariant subspace.
    /// \pre \ref compute() has been called.
    InnerQuantumState getBlockSize(BlockNumber in) const;

    /// Get the list of all Fock states spanning a given invariant subspace.
    /// \param[in] in Index of the invariant subspace.
    /// \pre \ref compute() has been called.
    std::vector<QuantumState> const& getFockStates(BlockNumber in) const;

    /// Get a specific Fock state from a given invariant subspace.
    /// \param[in] in Index of the invariant subspace.
    /// \param[in] i Index of the Fock state within the subspace.
    /// \pre \ref compute() has been called.
    QuantumState getFockState(BlockNumber in, InnerQuantumState i) const;

    /// Get the invariant subspace index a given Fock state belongs to.
    /// \param[in] in Fock state.
    /// \pre \ref compute() has been called.
    BlockNumber getBlockNumber(QuantumState in) const;

    /// For a given Fock state, get the index within the invariant subspace it belongs to.
    /// \param[in] in Fock state.
    /// \pre \ref compute() has been called.
    InnerQuantumState getInnerState(QuantumState in) const;

private:
    /// Initialize data members for a single un-partitioned Hilbert space.
    /// \param[in] Dim Dimension of the Hilbert space.
    void initSingleBlock(QuantumState Dim);
    /// Initialize data members for a partitioned Hilbert space.
    /// \param[in] partition Partition of the full Hilbert space into invariant subspaces.
    template <typename SpacePartitionType> void initMultipleBlocks(SpacePartitionType const& partition) {
        StateBlockIndex.resize(partition.dim());
        StatesContainer.resize(partition.n_subspaces(), std::vector<QuantumState>());
        foreach(partition, [this](QuantumState State, BlockNumber Block) {
            StateBlockIndex[State] = Block;
            StatesContainer[Block].push_back(State);
        })
            ;
    }
    /// Check if \ref compute() has already been called.
    void checkComputed() const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_STATESCLASSIFICATION_HPP
