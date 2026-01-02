//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/FieldOperatorContainer.hpp
/// \brief A container for creation and annihilation operators.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_HPP
#define POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_HPP

#include "Hamiltonian.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"

#include <set>
#include <unordered_map>

namespace Pomerol {

/// \addtogroup ED
///@{

/// \brief Container for instances of \ref CreationOperator and \ref  AnnihilationOperator.
///
/// This container class stores instances of \ref CreationOperator and \ref  AnnihilationOperator
/// in associative maps with keys being their respective single-particle indices.
/// It also provides methods that prepare and compute all stored \ref MonomialOperator objects
/// at once.
class FieldOperatorContainer {

    /// Storage of CreationOperator objects.
    std::unordered_map<ParticleIndex, CreationOperator> mapCreationOperators;
    /// Storage of AnnihilationOperator objects.
    std::unordered_map<ParticleIndex, AnnihilationOperator> mapAnnihilationOperators;

public:
    /// Constructor.
    /// \tparam IndexTypes Types of indices carried by the creation and annihilation operators.
    /// \param[in] IndexInfo Map for fermionic operator index tuples.
    /// \param[in] HS Hilbert space.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] in Set of all single-particle indices to store the \ref MonomialOperator objects
    ///            for. When empty, a set of all indices from IndexInfo is used.
    template <typename... IndexTypes>
    FieldOperatorContainer(IndexClassification<IndexTypes...> const& IndexInfo,
                           HilbertSpace<IndexTypes...> const& HS,
                           StatesClassification const& S,
                           Hamiltonian const& H,
                           std::set<ParticleIndex> in = {}) {
        if(in.empty()) {
            for(ParticleIndex p = 0; p < IndexInfo.getIndexSize(); ++p) {
                in.insert(p);
            }
        }
        for(auto p : in) {
            mapCreationOperators.emplace(p, CreationOperator(IndexInfo, HS, S, H, p));
            mapAnnihilationOperators.emplace(p, AnnihilationOperator(IndexInfo, HS, S, H, p));
        }
    }

    /// Prepare all stored creation and annihilation operators (allocate memory for them).
    /// \tparam IndexTypes Types of indices carried by the creation and annihilation operators.
    /// \param[in] HS Hilbert space.
    template <typename... IndexTypes> void prepareAll(HilbertSpace<IndexTypes...> const& HS) {
        for(auto& CX : mapCreationOperators)
            CX.second.prepare(HS);
        for(auto& C : mapAnnihilationOperators)
            C.second.prepare(HS);
    }

    /// Compute all stored creation and annihilation operators.
    /// \param[in] Tolerance Matrix elements with the absolute value equal or below this threshold
    ///                      are considered negligible.
    /// \pre \ref prepareAll() has been called.
    void computeAll(RealType Tolerance = 1e-8);

    /// Return a reference to a creation operator by its single-particle index.
    /// \param[in] in Single-particle index.
    CreationOperator const& getCreationOperator(ParticleIndex in) const;
    /// Return a reference to a annihilation operator by its single-particle index.
    /// \param[in] in Single-particle index.
    AnnihilationOperator const& getAnnihilationOperator(ParticleIndex in) const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_HPP
