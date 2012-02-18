//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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


#include "Vertex4.h"
#include <Eigen/LU> 

namespace Pomerol{

Vertex4Element::Vertex4Element(const TwoParticleGFContainer& Chi, const GFContainer& g,
                               const FourIndexObject::IndexCombination& Indices) :
    Thermal(Chi), Chi(Chi), g(g), Indices(Indices), Computed(false)
{}

void Vertex4Element::compute(long NumberOfMatsubaras)
{
    Storage.fill(this,NumberOfMatsubaras);
}

ComplexType Vertex4Element::value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    ComplexType Value = Chi(Indices,MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta*  g(Indices.Index1,Indices.Index3,MatsubaraNumber1)*
                        g(Indices.Index2,Indices.Index4,MatsubaraNumber2);
    if(MatsubaraNumber2 == MatsubaraNumber3)
        Value -= beta*  g(Indices.Index1,Indices.Index4,MatsubaraNumber1)*
                        g(Indices.Index2,Indices.Index3,MatsubaraNumber2);
    return Value;
}

ComplexType Vertex4Element::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    return Storage(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}

bool Vertex4Element::isComputed(void) const
{
    return Computed;
}

bool Vertex4Element::isVanishing(void) const
{
    // We need a smarter mechanism to detect this.
    return false;
}

Vertex4::Vertex4(const IndexClassification &IndexInfo, const TwoParticleGFContainer &Chi, const GFContainer &g) :
    FourIndexContainerObject(IndexInfo), Thermal(Chi), Chi(Chi), g(g)
{}

static Permutation4 EquivalentPermutations[4] =
    {permutations4[0],permutations4[1],permutations4[6],permutations4[7]};

void Vertex4::prepare(void)
{
    // TODO: This function is to be rewritten when the symmetry analyzer is done.
    ParticleIndex MaxIndex = IndexInfo.getIndexSize();
    for(ParticleIndex Index1=0; Index1 < MaxIndex; ++Index1)
    for(ParticleIndex Index2=Index1; Index2 < MaxIndex; ++Index2)
        for(ParticleIndex Index3=0; Index3 < MaxIndex; ++Index3)
        for(ParticleIndex Index4=Index3; Index4 < MaxIndex; ++Index4){
            Vertex4Element* tempV4E = new Vertex4Element(Chi,g,IndexCombination(Index1,Index2,Index3,Index4));
            ParticleIndex PI[4] = {Index1,Index2,Index3,Index4};
            for(size_t p=0; p<4; ++p){
                IndexCombination PermutedIndices(
                PI[EquivalentPermutations[p].perm[0]],
                PI[EquivalentPermutations[p].perm[1]],
                PI[EquivalentPermutations[p].perm[2]],
                PI[EquivalentPermutations[p].perm[3]]);
                ElementsMap[PermutedIndices] = new ElementWithPermFreq(tempV4E,EquivalentPermutations[p]);
                DEBUG("Adding Gamma4 component " << PermutedIndices <<
                       " with frequency permutation " << EquivalentPermutations[p] << ":" <<
                       " Vertex4Element at " << ElementsMap[PermutedIndices]->Element
                    )
            }
        }
}

void Vertex4::prepare(const std::set<IndexCombination>& InitialCombinations)
{
    // TODO: This function is to be rewritten when the symmetry analyzer is done.
    for (std::set<IndexCombination>::iterator it1=InitialCombinations.begin(); it1!=InitialCombinations.end(); ++it1){
        Vertex4Element* tempV4E = new Vertex4Element(Chi,g,*it1);
        ParticleIndex PI[4] = {(*it1).Index1,(*it1).Index2,(*it1).Index3,(*it1).Index4};
        for(size_t p=0; p<sizeof(EquivalentPermutations); ++p){
            IndexCombination PermutedIndices(
                PI[EquivalentPermutations[p].perm[0]],
                PI[EquivalentPermutations[p].perm[1]],
                PI[EquivalentPermutations[p].perm[2]],
                PI[EquivalentPermutations[p].perm[3]]);
            if(ElementsMap.count(PermutedIndices)==0){
                ElementsMap[PermutedIndices] = new ElementWithPermFreq(tempV4E,EquivalentPermutations[p]);
                DEBUG("Adding Gamma4 component " << PermutedIndices <<
                      " with frequency permutation " << EquivalentPermutations[p] << ":" <<
                      " Vertex4Element at " << ElementsMap[PermutedIndices]
                )
            }
        }
    }
}
 
void Vertex4::compute(long NumberOfMatsubaras)
{
#ifndef pomerolOpenMP
    for (std::map<IndexCombination,ElementWithPermFreq*>::iterator it1=ElementsMap.begin();it1!=ElementsMap.end();++it1){
        if(!it1->second->Element->isComputed())
            it1->second->Element->compute(NumberOfMatsubaras);
    };
#else
    std::vector<Vertex4Element*> items;
    for (std::map<IndexCombination,ElementWithPermFreq*>::iterator it1=ElementsMap.begin();it1!=ElementsMap.end();++it1){
        if (it1->second->Element->isComputed()) items.push_back(it1->second->Element);
    }
#pragma omp parallel for 
    for (int i = 0; i < (int) items.size(); ++i){
       if(!items[i]->isComputed())
            items[i]->compute(NumberOfMatsubaras);
    };
#endif
}

} // end of namespace Pomerol
