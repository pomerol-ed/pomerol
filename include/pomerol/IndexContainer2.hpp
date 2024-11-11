//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/IndexContainer2.hpp
/// \brief A CRTP base for container types whose elements are addressable by
/// two single-particle indices.
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_HPP
#define POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_HPP

#include "Index.hpp"
#include "IndexClassification.hpp"

#include <map>
#include <memory>
#include <set>
#include <utility>

namespace Pomerol {

/// \addtogroup Misc
///@{

/// \brief Base class for sparse container types whose elements are addressable
/// by two single-particle indices.
///
/// \tparam ElementType Type of an element.
/// \tparam SourceObject Type of the source object used to create the elements.
template <typename ElementType, typename SourceObject> class IndexContainer2 {
protected:
    /// Each of the two indices can change in the range [0; NumIndices[.
    ParticleIndex NumIndices;
    /// Sparse storage for the elements.
    std::map<IndexCombination2, std::shared_ptr<ElementType>> ElementsMap;

    /// Stored elements are created by calling Source.createElement(Indices).
    SourceObject const& Source;

    /// Generate a complete set of index combinations usable to address
    /// elements in the container.
    std::set<IndexCombination2> enumerateIndices() const;

public:
    /// Construct from a source object and an index classification object.
    /// The container is initially empty and shall be populated with elements
    /// by a subsequent call to \ref fill().
    /// \tparam IndexTypes Types of indices carried by a single creation/annihilation operator.
    /// \param[in] Source Source object used to create stored elements.
    /// \param[in] IndexInfo Classification of single-particle indices.
    template <typename... IndexTypes>
    IndexContainer2(SourceObject const& Source, IndexClassification<IndexTypes...> const& IndexInfo)
        : NumIndices(IndexInfo.getIndexSize()), Source(Source) {}

    /// Fill the container with elements from the source object.
    /// Each element is created by calling Source.createElement(IndexCombination).
    /// \param[in] Indices Set of index combinations of the elements to be created.
    ///            An empty set results in creation of elements for all possible index combinations.
    void fill(std::set<IndexCombination2> Indices = std::set<IndexCombination2>());

    /// Create a stored element from the source object by its index combination.
    /// \param[in] Indices Index combination of the element to be created.
    /// \return Reference to the created element.
    ElementType& create(IndexCombination2 const& Indices);

    /// Check if an element for a given index combination is stored in the container.
    /// \param[in] Indices Index combination.
    bool isInContainer(IndexCombination2 const& Indices) const;
    /// Check if an element for a given index combination is stored in the container.
    /// \param[in] Index1 First index in the combination.
    /// \param[in] Index2 Second index in the combination.
    bool isInContainer(ParticleIndex Index1, ParticleIndex Index2) const;

    /// Get a reference to a stored element by its index combination.
    /// \param[in] Indices Index combination.
    ElementType& operator()(IndexCombination2 const& Indices);
    /// Get a reference to a stored element by its index combination.
    /// \param[in] Index1 First index in the combination.
    /// \param[in] Index2 Second index in the combination.
    ElementType& operator()(ParticleIndex Index1, ParticleIndex Index2);
};

///@}

/////////////////////
// IndexContainer2 //
/////////////////////
template <typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType, SourceObject>::isInContainer(IndexCombination2 const& Indices) const {
    return ElementsMap.count(Indices) > 0;
}

template <typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType, SourceObject>::isInContainer(ParticleIndex Index1, ParticleIndex Index2) const {
    return isInContainer(IndexCombination2(Index1, Index2));
}

template <typename ElementType, typename SourceObject>
void IndexContainer2<ElementType, SourceObject>::fill(std::set<IndexCombination2> Indices) {
    // TODO: this method should use symmetry information
    // Indices should be split into equivalence classes with
    // the equivalence relation provided by a symmetry analyzer
    // The resulting classes may be extended afterwards (optionally)

    // remove existing elements
    ElementsMap.clear();

    std::set<IndexCombination2> II = Indices.empty() ? enumerateIndices() : std::move(Indices);

    for(auto const& ic : II) {
        if(!isInContainer(ic)) {
            create(ic);
        }
    }
}

template <typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType, SourceObject>::create(IndexCombination2 const& Indices) {
    std::shared_ptr<ElementType> pElement(Source.createElement(Indices));
    ElementsMap[Indices] = pElement;

    DEBUG("IndexContainer2::create() at " << this
                                          << ": "
                                             "added an element with indices "
                                          << Indices << " (" << pElement << ").");

    return *pElement;
}

template <typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType, SourceObject>::operator()(IndexCombination2 const& Indices) {
    auto iter = ElementsMap.find(Indices);

    if(iter == ElementsMap.end()) {
        DEBUG("IndexContainer2 at " << this << ": "
                                    << "cache miss for Index1=" << Indices.Index1 << ", Index2=" << Indices.Index2
                                    << "; add a new element to the container using source " << &Source);
        return create(Indices);
    }
    return *(iter->second);
}

template <typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType, SourceObject>::operator()(ParticleIndex Index1, ParticleIndex Index2) {
    return operator()(IndexCombination2(Index1, Index2));
}

template <typename ElementType, typename SourceObject>
inline std::set<IndexCombination2> IndexContainer2<ElementType, SourceObject>::enumerateIndices() const {
    std::set<IndexCombination2> AllIndices;

    for(ParticleIndex Index1 = 0; Index1 < NumIndices; ++Index1)
        for(ParticleIndex Index2 = 0; Index2 < NumIndices; ++Index2)
            AllIndices.emplace(Index1, Index2);

    return AllIndices;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_HPP
