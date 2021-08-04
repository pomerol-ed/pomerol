/** \file include/pomerol/IndexContainer2.h
** \brief A parent abstract class for 2-index objects
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_H
#define POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_H

#include "Index.hpp"
#include "IndexClassification.hpp"

#include <map>
#include <memory>
#include <set>
#include <utility>

namespace Pomerol {

/** A container to store components of a function with 2 indices. */
template<typename ElementType, typename SourceObject>
class IndexContainer2 {
protected:
    ParticleIndex NumIndices;
    std::map<IndexCombination2, std::shared_ptr<ElementType>> ElementsMap;

    SourceObject const& Source;

    std::set<IndexCombination2> enumerateInitialIndices() const;

public:

    template<typename... IndexTypes>
    IndexContainer2(SourceObject const& Source, const IndexClassification<IndexTypes...>& IndexInfo) :
        NumIndices(IndexInfo.getIndexSize()), Source(Source)
    {}

    void fill(std::set<IndexCombination2> InitialIndices = std::set<IndexCombination2>());
    ElementType& set(const IndexCombination2& Indices);

    bool isInContainer(const IndexCombination2& Indices) const;
    bool isInContainer(ParticleIndex Index1, ParticleIndex Index2) const;

    ElementType& operator()(const IndexCombination2& Indices);
    ElementType& operator()(ParticleIndex Index1, ParticleIndex Index2);
};

/////////////////////
// IndexContainer2 //
/////////////////////
template<typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType,SourceObject>::isInContainer(const IndexCombination2& Indices) const
{
    return ElementsMap.count(Indices) > 0;
}

template<typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType,SourceObject>::isInContainer(ParticleIndex Index1, ParticleIndex Index2) const
{
    return isInContainer(IndexCombination2(Index1, Index2));
}

template<typename ElementType, typename SourceObject>
void IndexContainer2<ElementType,SourceObject>::fill(std::set<IndexCombination2> InitialIndices)
{
    // TODO: this method should use symmetry information
    // InitialIndices should be split into equivalence classes with
    // the equivalence relation provided by a symmetry analyzer
    // The resulting classes may be extended afterwards (optionally)

    // remove existing elements
    ElementsMap.clear();

    std::set<IndexCombination2> II = InitialIndices.empty() ?
                                     enumerateInitialIndices() :
                                     std::move(InitialIndices);

    for(auto const& ic : II) {
        if(!isInContainer(ic)) {
            set(ic);
        }
    }
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType,SourceObject>::set(const IndexCombination2& Indices)
{
    std::shared_ptr<ElementType> pElement(Source.createElement(Indices));
    ElementsMap[Indices] = pElement;

    DEBUG("IndexContainer2::set() at " << this << ": "
          "added an element with indices " << Indices << " (" << pElement << ").");

    return *pElement;
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType,SourceObject>::operator()(const IndexCombination2& Indices)
{
    auto iter = ElementsMap.find(Indices);

    if(iter == ElementsMap.end()) {
        DEBUG("IndexContainer2 at " << this << ": " <<
              "cache miss for Index1=" << Indices.Index1 <<
              ", Index2=" << Indices.Index2 <<
              "; add a new element to the container using source " << &Source
        );
        return set(Indices);
    }
    return *(iter->second);
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType, SourceObject>::operator()
    (ParticleIndex Index1, ParticleIndex Index2)
{
    return operator()(IndexCombination2(Index1, Index2));
}

template<typename ElementType, typename SourceObject>
inline
std::set<IndexCombination2> IndexContainer2<ElementType,SourceObject>::enumerateInitialIndices() const
{
    std::set<IndexCombination2> AllIndices;

    for(ParticleIndex Index1 = 0; Index1 < NumIndices; ++Index1)
        for(ParticleIndex Index2 = 0; Index2 < NumIndices; ++Index2)
            AllIndices.emplace(Index1, Index2);

    return AllIndices;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEXCONTAINER2_H
