/** \file include/pomerol/IndexContainer2.h
** \brief A parent abstract class for 2-index objects
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_INDEXCONTAINER2_H
#define __INCLUDE_INDEXCONTAINER2_H

#include"IndexClassification.h"

#include<set>
#include<boost/shared_ptr.hpp>

namespace Pomerol {

/** A container to store components of a function with 2 indices. */
template<typename ElementType, typename SourceObject>
class IndexContainer2 {
protected:
    const IndexClassification& IndexInfo;
    std::map<IndexCombination2,boost::shared_ptr<ElementType> > ElementsMap;

    SourceObject* pSource;

    const std::set<IndexCombination2> enumerateInitialIndices(void) const;

public:

    IndexContainer2<ElementType,SourceObject>(SourceObject* pSource, const IndexClassification& IndexInfo);

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
IndexContainer2<ElementType,SourceObject>::IndexContainer2(SourceObject* pSource, const IndexClassification& IndexInfo) :
    IndexInfo(IndexInfo), pSource(pSource)
{}

template<typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType,SourceObject>::isInContainer(const IndexCombination2& Indices) const
{
    return ElementsMap.count(Indices) > 0;
}

template<typename ElementType, typename SourceObject>
bool IndexContainer2<ElementType,SourceObject>::isInContainer(ParticleIndex Index1, ParticleIndex Index2) const
{
    return isInContainer(IndexCombination2(Index1,Index2));
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

    std::set<IndexCombination2> II;
    if(InitialIndices.size()==0)           // If there are no indices provided,
        II = enumerateInitialIndices(); // Enumerate all possible combinations. 
    else
        II = InitialIndices;    // Otherwise use provided indices.

    for(typename std::set<IndexCombination2>::iterator iter = II.begin();
        iter != II.end(); iter++){
        if(!isInContainer(*iter)) {
            set(*iter);
            }
    }
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType,SourceObject>::set(const IndexCombination2& Indices)
{
    boost::shared_ptr<ElementType> pElement(pSource->createElement(Indices));
    ElementsMap[Indices] = pElement;

    DEBUG("IndexContainer2::set() at " << this << ": "
          "added an element with indices " << Indices << " (" << pElement << ").");

    return *pElement;
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType,SourceObject>::operator()(const IndexCombination2& Indices)
{
    typename std::map<IndexCombination2,boost::shared_ptr<ElementType> >::iterator
        iter = ElementsMap.find(Indices);

    if(iter == ElementsMap.end()){
        DEBUG("IndexContainer2 at " << this << ": " <<
              "cache miss for Index1=" << Indices.Index1 <<
              ", Index2=" << Indices.Index2 <<
              "; add a new element to the container using source " << pSource
        );
        return set(Indices);
    }
    return *(iter->second);
}

template<typename ElementType, typename SourceObject>
ElementType& IndexContainer2<ElementType,SourceObject>::operator()
    (ParticleIndex Index1, ParticleIndex Index2)
{
    return operator()(IndexCombination2(Index1,Index2));
}

template<typename ElementType, typename SourceObject>
inline
const std::set<IndexCombination2> IndexContainer2<ElementType,SourceObject>::enumerateInitialIndices(void) const
{
    std::set<IndexCombination2> AllIndices;

    ParticleIndex Size = IndexInfo.getIndexSize();
    for(ParticleIndex Index1=0; Index1<Size; ++Index1)
    for(ParticleIndex Index2=0; Index2<Size; ++Index2)
        AllIndices.insert(IndexCombination2(Index1,Index2));

    return AllIndices;
}


} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_INDEXCONTAINER2_H
