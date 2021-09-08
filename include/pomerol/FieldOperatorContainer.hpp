//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/FieldOperatorContainer.h
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_H
#define POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_H

#include "Hamiltonian.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"

#include <set>
#include <unordered_map>

namespace Pomerol {

/** This class represents a container to store and retrieve FieldOperators ( CreationOperator or AnnihilationOperator
 * rotated to eigenvector basis of Hamiltonian H ) for a given Index.
 * If no field operator is yet initialized then calculation of the field operator is done.
 */
class FieldOperatorContainer {

    /** A map which gives a link to the CreationOperator for a given index */
    std::unordered_map<ParticleIndex, CreationOperator> mapCreationOperators;
    /** A map which gives a link to the AnnihilationOperator for a given index */
    std::unordered_map<ParticleIndex, AnnihilationOperator> mapAnnihilationOperators;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] IndexInfo A reference to a IndexClassification
     */
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

    template <typename... IndexTypes> void prepareAll(HilbertSpace<IndexTypes...> const& HS) {
        for(auto& CX : mapCreationOperators)
            CX.second.prepare(HS);
        for(auto& C : mapAnnihilationOperators)
            C.second.prepare(HS);
    }
    void computeAll();

    /** Returns the CreationOperator for a given Index. Makes on-demand computation. */
    CreationOperator const& getCreationOperator(ParticleIndex in) const;
    /** Returns the AnnihilationOperator for a given Index. Makes on-demand computation */
    AnnihilationOperator const& getAnnihilationOperator(ParticleIndex in) const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_FIELDOPERATORCONTAINER_H
