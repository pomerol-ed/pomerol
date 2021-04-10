 //
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


/** \file src/IndexContainer4.h
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_INDEXCONTAINER4_H
#define __INCLUDE_INDEXCONTAINER4_H

#include"IndexClassification.h"

#include <memory>
#include <set>

namespace Pomerol {

/** Decorator around an element to intercept operator(n1,n2,n3) calls nd permute Matsubara frequencies. */
template<typename ElementType>
struct ElementWithPermFreq
{
    std::shared_ptr<ElementType> pElement;
    const Permutation4 FrequenciesPermutation;

    ElementWithPermFreq(std::shared_ptr<ElementType> pElement, Permutation4 FrequenciesPermutation);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    operator ElementType&();
};

/** A container to store components of a function with 4 indices. */
template<typename ElementType, typename SourceObject>
class IndexContainer4 {
protected:
    ParticleIndex NumIndices;

    SourceObject* pSource;

    const std::set<IndexCombination4> enumerateInitialIndices(void) const;

public:
    std::map<IndexCombination4,ElementWithPermFreq<ElementType> > ElementsMap;
    std::map<IndexCombination4, std::shared_ptr<ElementType>  > NonTrivialElements;

    template<typename... IndexTypes>
    IndexContainer4(SourceObject* pSource, const IndexClassification<IndexTypes...>& IndexInfo) :
        NumIndices(IndexInfo.getIndexSize()), pSource(pSource)
    {}

    void fill(std::set<IndexCombination4> InitialIndices = std::set<IndexCombination4>());
    ElementWithPermFreq<ElementType>& set(const IndexCombination4& Indices);

    bool isInContainer(const IndexCombination4& Indices) const;
    bool isInContainer( ParticleIndex Index1, ParticleIndex Index2,
                        ParticleIndex Index3, ParticleIndex Index4) const;

    ElementWithPermFreq<ElementType>& operator()(const IndexCombination4& Indices);
    ElementWithPermFreq<ElementType>& operator()(ParticleIndex Index1, ParticleIndex Index2,
                                                ParticleIndex Index3, ParticleIndex Index4);
};

/////////////////////////
// ElementWithPermFreq //
/////////////////////////
template<typename ElementType> inline
ElementWithPermFreq<ElementType>::ElementWithPermFreq(
    std::shared_ptr<ElementType> pElement, Permutation4 FrequenciesPermutation) :
    pElement(pElement), FrequenciesPermutation(FrequenciesPermutation)
{}

template<typename ElementType> inline
ComplexType ElementWithPermFreq<ElementType>::operator()(
    long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long MatsubaraNumbers[4] = {MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3,
                                MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3};
    return (*pElement)(MatsubaraNumbers[FrequenciesPermutation.perm[0]],
                       MatsubaraNumbers[FrequenciesPermutation.perm[1]],
                       MatsubaraNumbers[FrequenciesPermutation.perm[2]])
                    *RealType(FrequenciesPermutation.sign);
}

template<typename ElementType> inline
ElementWithPermFreq<ElementType>::operator ElementType&()
{
    return *pElement;
}

/////////////////////
// IndexContainer4 //
/////////////////////
template<typename ElementType, typename SourceObject>
inline
bool IndexContainer4<ElementType,SourceObject>::isInContainer(const IndexCombination4& Indices) const
{
    return ElementsMap.count(Indices) > 0;
}

template<typename ElementType, typename SourceObject>
inline
bool IndexContainer4<ElementType,SourceObject>::isInContainer(  ParticleIndex Index1, ParticleIndex Index2,
                                                                ParticleIndex Index3, ParticleIndex Index4) const
{
    return isInContainer(IndexCombination4(Index1,Index2,Index3,Index4));
}

template<typename ElementType, typename SourceObject>
inline
void IndexContainer4<ElementType,SourceObject>::fill(std::set<IndexCombination4> InitialIndices)
{
    // TODO: this method should use symmetry information
    // InitialIndices should be split into equivalence classes with
    // the equivalence relation provided by a symmetry analyzer
    // The resulting classes may be extended afterwards (optionally)

    // remove existing elements
    ElementsMap.clear();

    std::set<IndexCombination4> II;
    if(InitialIndices.size()==0)           // If there are no indices provided,
        II = enumerateInitialIndices(); // Enumerate all possible combinations.
    else
        II = InitialIndices;    // Otherwise use provided indices.

    for(typename std::set<IndexCombination4>::iterator iter = II.begin();
        iter != II.end(); iter++){
        if(!isInContainer(*iter))
            set(*iter);
    }
}

template<typename ElementType, typename SourceObject>
inline
ElementWithPermFreq<ElementType>& IndexContainer4<ElementType,SourceObject>::set(const IndexCombination4& Indices)
{
    // TODO: rewrite this method entirely (merge with fill()?)
    std::shared_ptr<ElementType> pElement(pSource->createElement(Indices));
    typename std::map<IndexCombination4,ElementWithPermFreq<ElementType> >::iterator iter =
        ElementsMap.insert(
            std::pair<IndexCombination4,ElementWithPermFreq<ElementType> >
                (Indices,ElementWithPermFreq<ElementType>(pElement,permutations4[0]))).first;

    DEBUG("IndexContainer4::fill() at " << this << ": " <<
        "added an element with indices " << Indices <<
        " and frequency permutation " << permutations4[0] <<
        " (" << pElement << ").");

    bool SameCIndices = (Indices.Index1==Indices.Index2);
    bool SameCXIndices = (Indices.Index3==Indices.Index4);

    NonTrivialElements.insert(std::make_pair(Indices,pElement));

    if(!SameCIndices){
        IndexCombination4 Indices2134(Indices.Index2,Indices.Index1,Indices.Index3,Indices.Index4);
        if(!isInContainer(Indices2134)){
            ElementsMap.insert(
            std::pair<IndexCombination4,ElementWithPermFreq<ElementType> >(
                Indices2134,
                ElementWithPermFreq<ElementType>(pElement,permutations4[6])));
            DEBUG("IndexContainer4::fill() at " << this << ": " <<
                "added an element with indices " << Indices <<
                " and frequency permutation " << permutations4[6] <<
                " (" << pElement << ").");
        }
    }
    if(!SameCXIndices){
        IndexCombination4 Indices1243(Indices.Index1,Indices.Index2,Indices.Index4,Indices.Index3);
        if(!isInContainer(Indices1243)){
            ElementsMap.insert(
                std::pair<IndexCombination4,ElementWithPermFreq<ElementType> >(
                Indices1243,
                ElementWithPermFreq<ElementType>(pElement,permutations4[1])));
            DEBUG("IndexContainer4::fill() at " << this << ": " <<
                "added an element with indices " << Indices <<
                " and frequency permutation " << permutations4[1] <<
                " (" << pElement << ").");
        }
    }
    if(!SameCIndices && !SameCXIndices){
        IndexCombination4 Indices2143(Indices.Index2,Indices.Index1,Indices.Index4,Indices.Index3);
        if(!isInContainer(Indices2143)){
            ElementsMap.insert(
                std::pair<IndexCombination4,ElementWithPermFreq<ElementType> >(
                Indices2143,
                ElementWithPermFreq<ElementType>(pElement,permutations4[7])));
            DEBUG("IndexContainer4::fill() at " << this << ": " <<
                "added an element with indices " << Indices <<
                " and frequency permutation " << permutations4[7] <<
                " (" << pElement << ").");
        }
    }

    return iter->second;
}

template<typename ElementType, typename SourceObject>
inline
ElementWithPermFreq<ElementType>& IndexContainer4<ElementType,SourceObject>::operator()
    (const IndexCombination4& Indices)
{
    typename std::map<IndexCombination4,ElementWithPermFreq<ElementType> >::iterator iter
        = ElementsMap.find(Indices);

    if(iter == ElementsMap.end()){
        DEBUG("IndexContainer4 at " << this << ": " <<
              "cache miss for Index1=" << Indices.Index1 <<
              ", Index2=" << Indices.Index2 <<
              ", Index3=" << Indices.Index3 <<
              ", Index4=" << Indices.Index4 <<
              "; add a new element to the container using source " << pSource
        );
        return set(Indices);
    }else
        return iter->second;
}

template<typename ElementType, typename SourceObject>
inline
ElementWithPermFreq<ElementType>& IndexContainer4<ElementType,SourceObject>::operator()
    (ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4)
{
    return operator()(IndexCombination4(Index1,Index2,Index3,Index4));
}

template<typename ElementType, typename SourceObject>
inline
const std::set<IndexCombination4> IndexContainer4<ElementType,SourceObject>::enumerateInitialIndices(void) const
{
    std::set<IndexCombination4> AllIndices;

    for(ParticleIndex Index1=0; Index1<NumIndices; ++Index1)
    for(ParticleIndex Index2=Index1; Index2<NumIndices; ++Index2)
        for(ParticleIndex Index3=0; Index3<NumIndices; ++Index3)
        for(ParticleIndex Index4=Index3; Index4<NumIndices; ++Index4)
            AllIndices.insert(IndexCombination4(Index1,Index2,Index3,Index4));

    return AllIndices;
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_INDEXCONTAINER4_H
