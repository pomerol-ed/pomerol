//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file src/IndexContainer4.h
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER4_HPP
#define POMEROL_INCLUDE_POMEROL_INDEXCONTAINER4_HPP

#include "Index.hpp"
#include "IndexClassification.hpp"

#include <map>
#include <memory>
#include <set>
#include <utility>

namespace Pomerol {

/** Decorator around an element to intercept operator(n1,n2,n3) calls nd permute Matsubara frequencies. */
template <typename ElementType> struct ElementWithPermFreq {
    std::shared_ptr<ElementType> pElement;
    Permutation4 const FrequenciesPermutation;

    ElementWithPermFreq(std::shared_ptr<ElementType> pElement, Permutation4 const& FrequenciesPermutation);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    operator ElementType&();
};

/** A container to store components of a function with 4 indices. */
template <typename ElementType, typename SourceObject> class IndexContainer4 {
protected:
    ParticleIndex NumIndices;

    SourceObject const& Source;

    std::set<IndexCombination4> enumerateInitialIndices() const;

public:
    std::map<IndexCombination4, ElementWithPermFreq<ElementType>> ElementsMap;
    std::map<IndexCombination4, std::shared_ptr<ElementType>> NonTrivialElements;

    template <typename... IndexTypes>
    IndexContainer4(SourceObject const& Source, IndexClassification<IndexTypes...> const& IndexInfo)
        : NumIndices(IndexInfo.getIndexSize()), Source(Source) {}

    void fill(std::set<IndexCombination4> InitialIndices = std::set<IndexCombination4>());
    ElementWithPermFreq<ElementType>& set(IndexCombination4 const& Indices);

    bool isInContainer(IndexCombination4 const& Indices) const;
    bool isInContainer(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4) const;

    ElementWithPermFreq<ElementType>& operator()(IndexCombination4 const& Indices);
    ElementWithPermFreq<ElementType>&
    operator()(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4);
};

/////////////////////////
// ElementWithPermFreq //
/////////////////////////
template <typename ElementType>
inline ElementWithPermFreq<ElementType>::ElementWithPermFreq(std::shared_ptr<ElementType> pElement,
                                                             Permutation4 const& FrequenciesPermutation)
    : pElement(pElement), FrequenciesPermutation(FrequenciesPermutation) {}

template <typename ElementType>
inline ComplexType ElementWithPermFreq<ElementType>::operator()(long MatsubaraNumber1,
                                                                long MatsubaraNumber2,
                                                                long MatsubaraNumber3) const {
    long MatsubaraNumbers[4] = {MatsubaraNumber1,
                                MatsubaraNumber2,
                                MatsubaraNumber3,
                                MatsubaraNumber1 + MatsubaraNumber2 - MatsubaraNumber3};
    return (*pElement)(MatsubaraNumbers[FrequenciesPermutation.perm[0]],
                       MatsubaraNumbers[FrequenciesPermutation.perm[1]],
                       MatsubaraNumbers[FrequenciesPermutation.perm[2]]) *
           RealType(FrequenciesPermutation.sign);
}

template <typename ElementType> inline ElementWithPermFreq<ElementType>::operator ElementType&() {
    return *pElement;
}

/////////////////////
// IndexContainer4 //
/////////////////////
template <typename ElementType, typename SourceObject>
inline bool IndexContainer4<ElementType, SourceObject>::isInContainer(IndexCombination4 const& Indices) const {
    return ElementsMap.count(Indices) > 0;
}

template <typename ElementType, typename SourceObject>
inline bool IndexContainer4<ElementType, SourceObject>::isInContainer(ParticleIndex Index1,
                                                                      ParticleIndex Index2,
                                                                      ParticleIndex Index3,
                                                                      ParticleIndex Index4) const {
    return isInContainer(IndexCombination4(Index1, Index2, Index3, Index4));
}

template <typename ElementType, typename SourceObject>
inline void IndexContainer4<ElementType, SourceObject>::fill(std::set<IndexCombination4> InitialIndices) {
    // TODO: this method should use symmetry information
    // InitialIndices should be split into equivalence classes with
    // the equivalence relation provided by a symmetry analyzer
    // The resulting classes may be extended afterwards (optionally)

    // remove existing elements
    ElementsMap.clear();

    std::set<IndexCombination4> II = InitialIndices.empty() ? enumerateInitialIndices() : std::move(InitialIndices);

    for(auto const& ic : II) {
        if(!isInContainer(ic)) {
            set(ic);
        }
    }
}

template <typename ElementType, typename SourceObject>
inline ElementWithPermFreq<ElementType>&
IndexContainer4<ElementType, SourceObject>::set(IndexCombination4 const& Indices) {
    std::shared_ptr<ElementType> pElement(Source.createElement(Indices));
    auto iter = ElementsMap.emplace(Indices, ElementWithPermFreq<ElementType>(pElement, permutations4[0])).first;

    DEBUG("IndexContainer4::fill() at " << this << ": "
                                        << "added an element with indices " << Indices << " and frequency permutation "
                                        << permutations4[0] << " (" << pElement << ").");

    bool SameCIndices = (Indices.Index1 == Indices.Index2);
    bool SameCXIndices = (Indices.Index3 == Indices.Index4);

    NonTrivialElements.emplace(Indices, pElement);

    if(!SameCIndices) {
        IndexCombination4 Indices2134(Indices.Index2, Indices.Index1, Indices.Index3, Indices.Index4);
        if(!isInContainer(Indices2134)) {
            ElementsMap.emplace(Indices2134, ElementWithPermFreq<ElementType>(pElement, permutations4[6]));
            DEBUG("IndexContainer4::fill() at " << this << ": "
                                                << "added an element with indices " << Indices
                                                << " and frequency permutation " << permutations4[6] << " (" << pElement
                                                << ").");
        }
    }
    if(!SameCXIndices) {
        IndexCombination4 Indices1243(Indices.Index1, Indices.Index2, Indices.Index4, Indices.Index3);
        if(!isInContainer(Indices1243)) {
            ElementsMap.emplace(Indices1243, ElementWithPermFreq<ElementType>(pElement, permutations4[1]));
            DEBUG("IndexContainer4::fill() at " << this << ": "
                                                << "added an element with indices " << Indices
                                                << " and frequency permutation " << permutations4[1] << " (" << pElement
                                                << ").");
        }
    }
    if(!SameCIndices && !SameCXIndices) {
        IndexCombination4 Indices2143(Indices.Index2, Indices.Index1, Indices.Index4, Indices.Index3);
        if(!isInContainer(Indices2143)) {
            ElementsMap.emplace(Indices2143, ElementWithPermFreq<ElementType>(pElement, permutations4[7]));
            DEBUG("IndexContainer4::fill() at " << this << ": "
                                                << "added an element with indices " << Indices
                                                << " and frequency permutation " << permutations4[7] << " (" << pElement
                                                << ").");
        }
    }

    return iter->second;
}

template <typename ElementType, typename SourceObject>
inline ElementWithPermFreq<ElementType>&
IndexContainer4<ElementType, SourceObject>::operator()(IndexCombination4 const& Indices) {
    auto iter = ElementsMap.find(Indices);

    if(iter == ElementsMap.end()) {
        DEBUG("IndexContainer4 at " << this << ": "
                                    << "cache miss for Index1=" << Indices.Index1 << ", Index2=" << Indices.Index2
                                    << ", Index3=" << Indices.Index3 << ", Index4=" << Indices.Index4
                                    << "; add a new element to the container using source " << &Source);
        return set(Indices);
    } else
        return iter->second;
}

template <typename ElementType, typename SourceObject>
inline ElementWithPermFreq<ElementType>& IndexContainer4<ElementType, SourceObject>::operator()(ParticleIndex Index1,
                                                                                                ParticleIndex Index2,
                                                                                                ParticleIndex Index3,
                                                                                                ParticleIndex Index4) {
    return operator()(IndexCombination4(Index1, Index2, Index3, Index4));
}

template <typename ElementType, typename SourceObject>
inline std::set<IndexCombination4> IndexContainer4<ElementType, SourceObject>::enumerateInitialIndices() const {
    std::set<IndexCombination4> AllIndices;

    for(ParticleIndex Index1 = 0; Index1 < NumIndices; ++Index1)
        for(ParticleIndex Index2 = Index1; Index2 < NumIndices; ++Index2)
            for(ParticleIndex Index3 = 0; Index3 < NumIndices; ++Index3)
                for(ParticleIndex Index4 = Index3; Index4 < NumIndices; ++Index4)
                    AllIndices.emplace(Index1, Index2, Index3, Index4);

    return AllIndices;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER4_HPP
