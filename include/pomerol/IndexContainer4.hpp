//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/IndexContainer4.hpp
/// \brief A CRTP base for container types whose elements are addressable by
/// four single-particle indices.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER4_HPP
#define POMEROL_INCLUDE_POMEROL_INDEXCONTAINER4_HPP

#include "Index.hpp"
#include "IndexClassification.hpp"

#include <map>
#include <memory>
#include <set>
#include <utility>

namespace Pomerol {

/// \addtogroup Misc
///@{

/// \brief A decorator that permutes indices of Matsubara frequencies
/// in calls to operator().
///
/// A decorator that intercepts calls to operator()(n1, n2, n3) and forwards them to an underlying object
/// while permuting the Matsubara frequency indices n1, n2, n3.
/// \tparam ElementType Underlying object type.
template <typename ElementType> struct ElementWithPermFreq {
    /// The underlying callable object.
    std::shared_ptr<ElementType> pElement;
    /// The permutation of the frequency indices.
    Permutation4 const FrequenciesPermutation;

    /// Constructor.
    /// \param[in] pElement The object to be decorated.
    /// \param[in] FrequenciesPermutation The Matsubara frequency index permutation.
    ElementWithPermFreq(std::shared_ptr<ElementType> pElement, Permutation4 const& FrequenciesPermutation);

    /// Call the underlying object with the permuted Matsubara frequency indices.
    /// \param[in] MatsubaraNumber1 First Matsubara frequency index.
    /// \param[in] MatsubaraNumber2 Second Matsubara frequency index.
    /// \param[in] MatsubaraNumber3 Third Matsubara frequency index.
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /// Return a reference to the underlying object.
    operator ElementType&();
};

/// \brief Base class for sparse container types whose elements are addressable
/// by four single-particle indices.
///
/// Base class for sparse container types whose elements are addressable by
/// four single-particle indices. The stored elements are also decorated by \ref ElementWithPermFreq.
/// \tparam ElementType Type of an element.
/// \tparam SourceObject Type of the source object used to create the elements.
template <typename ElementType, typename SourceObject> class IndexContainer4 {
protected:
    /// Each of the four indices can change in the range [0; NumIndices[.
    ParticleIndex NumIndices;

    /// Stored elements are created by calling Source.createElement(Indices).
    SourceObject const& Source;

    /// Generate a complete set of index combinations usable to address
    /// elements in the container.
    std::set<IndexCombination4> enumerateIndices() const;

    /// Sparse storage for the decorated elements.
    std::map<IndexCombination4, ElementWithPermFreq<ElementType>> ElementsMap;
    /// Sparse storage for the plain (non-decorated) elements.
    std::map<IndexCombination4, std::shared_ptr<ElementType>> NonTrivialElements;

public:
    /// Construct from a source object and an index classification object.
    /// The container is initially empty and shall be populated with elements
    /// by a subsequent call to \ref fill().
    /// \tparam IndexTypes Types of indices carried by a single creation/annihilation operator.
    /// \param[in] Source Source object used to create stored elements.
    /// \param[in] IndexInfo Classification of single-particle indices.
    template <typename... IndexTypes>
    IndexContainer4(SourceObject const& Source, IndexClassification<IndexTypes...> const& IndexInfo)
        : NumIndices(IndexInfo.getIndexSize()), Source(Source) {}

    /// Fill the container with elements from the source object.
    /// Each element is created by calling Source.createElement(IndexCombination).
    /// \param[in] Indices Set of index combinations of the elements to be created.
    ///            An empty set results in creation of elements for all possible index combinations.
    void fill(std::set<IndexCombination4> Indices = std::set<IndexCombination4>());

    /// Create a stored element from the source object by its index combination.
    /// \param[in] Indices Index combination of the element to be created.
    /// \return Reference to the created element.
    ElementWithPermFreq<ElementType>& create(IndexCombination4 const& Indices);

    /// Check if an element for a given index combination is stored in the container.
    /// \param[in] Indices Index combination.
    bool isInContainer(IndexCombination4 const& Indices) const;
    /// Check if an element for a given index combination is stored in the container.
    /// \param[in] Index1 First index in the combination.
    /// \param[in] Index2 Second index in the combination.
    /// \param[in] Index3 Third index in the combination.
    /// \param[in] Index4 Fourth index in the combination.
    bool isInContainer(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4) const;

    /// Get a reference to a stored element by its index combination.
    /// \param[in] Indices Index combination.
    ElementWithPermFreq<ElementType>& operator()(IndexCombination4 const& Indices);
    /// Get a reference to a stored element by its index combination.
    /// \param[in] Index1 First index in the combination.
    /// \param[in] Index2 Second index in the combination.
    /// \param[in] Index3 Third index in the combination.
    /// \param[in] Index4 Fourth index in the combination.
    ElementWithPermFreq<ElementType>&
    operator()(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4);
};

///@}

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
inline void IndexContainer4<ElementType, SourceObject>::fill(std::set<IndexCombination4> Indices) {
    // TODO: this method should use symmetry information
    // Indices should be split into equivalence classes with
    // the equivalence relation provided by a symmetry analyzer
    // The resulting classes may be extended afterwards (optionally)

    // remove existing elements
    ElementsMap.clear();

    std::set<IndexCombination4> II = Indices.empty() ? enumerateIndices() : std::move(Indices);

    for(auto const& ic : II) {
        if(!isInContainer(ic)) {
            create(ic);
        }
    }
}

template <typename ElementType, typename SourceObject>
inline ElementWithPermFreq<ElementType>&
IndexContainer4<ElementType, SourceObject>::create(IndexCombination4 const& Indices) {
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
        return create(Indices);
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
inline std::set<IndexCombination4> IndexContainer4<ElementType, SourceObject>::enumerateIndices() const {
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
