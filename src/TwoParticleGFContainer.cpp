//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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


#include "TwoParticleGFContainer.h"

namespace Pomerol{

extern std::ostream& OUTPUT_STREAM;
/*=========================================================================*/

TwoParticleGFContainer::Element::Element(TwoParticleGF* Computable2PGF, Permutation4 FrequenciesPermutation, bool isComputed):
    Computable2PGF(Computable2PGF),FrequenciesPermutation(FrequenciesPermutation),isComputed(isComputed)
{}

FourIndexContainerObject::IndexCombination TwoParticleGFContainer::Element::getLinkedIndices()
{
    return FourIndexContainerObject::IndexCombination(Computable2PGF->getIndex(0),Computable2PGF->getIndex(1),Computable2PGF->getIndex(2),Computable2PGF->getIndex(3));
}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,IndexClassification& IndexInfo, FieldOperatorContainer& Operators):
    ComputableObject(),Thermal(DM),S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators),NumberOfMatsubaras(0)
{
};

void TwoParticleGFContainer::readInitialIndices(std::vector<IndexCombination*>& in)
{
InitialCombinations=in;
for (std::vector<IndexCombination*>::const_iterator it1=in.begin(); it1!=in.end(); ++it1){
}
};

void TwoParticleGFContainer::defineInitialIndices()
{
if (InitialCombinations.size()==0)
for (ParticleIndex i1=0; i1<IndexInfo.getIndexSize(); ++i1)
    for (ParticleIndex i2=0; i2<IndexInfo.getIndexSize(); ++i2)
        for (ParticleIndex i3=0; i3<IndexInfo.getIndexSize(); ++i3){
            ParticleIndex i4=i1+i2;
            if (i3<=i1+i2){ // ParticleIndex is an unsigned short, so no negative value should be put in
                i4=i1+i2-i3;
                if (i4<IndexInfo.getIndexSize()){
                    IndexCombination *current = new IndexCombination(i1,i2,i3,i4);
                        InitialCombinations.push_back(current);
                    };
                };
            };
};

void TwoParticleGFContainer::prepare()
{
if (Status < Prepared){
    defineInitialIndices();
    for (std::vector<IndexCombination*>::const_iterator it1=InitialCombinations.begin(); it1!=InitialCombinations.end(); ++it1) if (mapNonTrivialCombinations.count(**it1)==0 ){
        AnnihilationOperator &C1 = Operators.getAnnihilationOperator((*it1)->Indices[0]);
        AnnihilationOperator &C2 = Operators.getAnnihilationOperator((*it1)->Indices[1]);
        CreationOperator     &CX3 = Operators.getCreationOperator   ((*it1)->Indices[2]);
        CreationOperator     &CX4 = Operators.getCreationOperator   ((*it1)->Indices[3]);
        TwoParticleGF * temp2PGF = new TwoParticleGF(S,H,C1,C2,CX3,CX4,DM);
        temp2PGF->prepare();
        if (!temp2PGF->vanishes()){
            Element *ComputableElement = new Element(temp2PGF,TrivialOperatorPermutations[0],true); 
            mapNonTrivialCombinations[**it1] = ComputableElement;
            NonTrivialCombinations.push_back(*it1);
            INFO("TwoParticleGFContainer: assigned " << **it1 << " TwoParticleGF for calculation");
            for (std::list<IndexClassification::IndexPermutation>::const_iterator index_perm = IndexInfo.getEquivalentIndexPermutations().begin();
                index_perm!=IndexInfo.getEquivalentIndexPermutations().end(); ++index_perm){
                    IndexCombination* permuted = new IndexCombination(index_perm->Permutation[(*it1)->Indices[0]],
                                                                      index_perm->Permutation[(*it1)->Indices[1]],
                                                                      index_perm->Permutation[(*it1)->Indices[2]],
                                                                      index_perm->Permutation[(*it1)->Indices[3]]
                                                                     );
                    Element* permutedElement = new Element(temp2PGF,TrivialOperatorPermutations[0], false);
                    if (mapNonTrivialCombinations.count(*permuted)==0 || 
                    (mapNonTrivialCombinations.count(*permuted)>0 && mapNonTrivialCombinations[*permuted]->isComputed && *permuted!= mapNonTrivialCombinations[*permuted]->getLinkedIndices())){

                        mapNonTrivialCombinations[*permuted] = permutedElement;
                        NonTrivialCombinations.push_back(permuted);
                        INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << "->" << **it1 );
                        INFO(" with a frequencies permutation 1234= 1234 due to symmetry considerations");
                        addInnerPermutationsOfIndexCombination(permuted);
                        }
                    else{
                        INFO("TwoParticleGFContainer: Detected frequency symmetry in " << *permuted << " : 1234=" << mapNonTrivialCombinations[*permuted]->FrequenciesPermutation << "!")
                        };
                };
            addInnerPermutationsOfIndexCombination(*it1);
                
        };
    };
    Status = Prepared;
    InitialCombinations.clear();
    INFO("");
    INFO("TwoParticleGCContainer: Following combinations will be obtained:")
    for (std::vector<IndexCombination*>::const_iterator it1=NonTrivialCombinations.begin(); it1!=NonTrivialCombinations.end(); ++it1){
        INFO_NONEWLINE(**it1 << " calculated: ");
        if (mapNonTrivialCombinations[**it1]->isComputed) INFO("y")
        else {INFO("n, is linked to " << mapNonTrivialCombinations[**it1]->getLinkedIndices()); 
             };
        }
    INFO("");
};
};

void TwoParticleGFContainer::addInnerPermutationsOfIndexCombination(const IndexCombination *in)
{
    IndexCombination in_computed = IndexCombination(mapNonTrivialCombinations[*in]->getLinkedIndices());
    for (int perm = 1; perm < 4; ++perm){
        Permutation4 perm4 = TrivialOperatorPermutations[perm];
        IndexCombination *permuted = new IndexCombination(in->Indices[perm4.perm[0]],in->Indices[perm4.perm[1]],in->Indices[perm4.perm[2]],in->Indices[perm4.perm[3]]);
        if (*permuted!=*in){
            if (mapNonTrivialCombinations.count(*permuted)==0 || 
                (mapNonTrivialCombinations.count(*permuted)>0 && mapNonTrivialCombinations[*permuted]->isComputed && *permuted!= mapNonTrivialCombinations[*permuted]->getLinkedIndices())){
                Element* permutedElement = new Element(mapNonTrivialCombinations[in_computed]->Computable2PGF,perm4,false);
                mapNonTrivialCombinations[*permuted] = permutedElement;
                NonTrivialCombinations.push_back(permuted);
                INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << "->" << in_computed);
                INFO_NONEWLINE(" with a frequencies permutation 1234=" << perm4);
                INFO(" due to internal indices permutation");
            }
            else if (mapNonTrivialCombinations[*permuted]->FrequenciesPermutation != TrivialOperatorPermutations[0]){
                // A thin moment. One checks that current permutation was not entered before i.e. by symmetry considerations or as an initial combination.
                // It might happen, that an obtained permutation is an exact copy not of IndexCombination *in, but rather the one
                // in NonTrivialCombinations, entered before. This 'if' operator is introduced to avoid stupid messages like 
                // "TwoParticleGFContainer: Detected frequency symmetry in ... : 1234=1234"
                INFO_NONEWLINE("TwoParticleGFContainer: Detected frequency symmetry in " << *permuted << " : 1234=");
                INFO(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation << "!");
                };

        };
    };
}

void TwoParticleGFContainer::compute(long NumberOfMatsubaras)
{
if (Status == Prepared){
    this->NumberOfMatsubaras = NumberOfMatsubaras;
#ifndef pomerolOpenMP
    for (std::map<IndexCombination,Element*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
        it1->second->Computable2PGF->compute(NumberOfMatsubaras);
        };
#else
    std::vector<TwoParticleGF*> items;
    for (std::map<IndexCombination,Element*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
        if (it1->second->isComputed) items.push_back(it1->second->Computable2PGF);
        }
    #pragma omp parallel for 
    for (int i = 0; i < (int) items.size(); ++i){
        items[i]->compute(NumberOfMatsubaras);
        };
#endif
    Status = Computed;    
    };
};


ComplexType TwoParticleGFContainer::operator()(const IndexCombination& in, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    if ( mapNonTrivialCombinations.count(in) > 0 ){
        Element* where = mapNonTrivialCombinations[in];
        long MatsubaraNumbers[4] = {MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3,MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3};
        ComplexType Value = (*where->Computable2PGF)(
                                                    MatsubaraNumbers[where->FrequenciesPermutation.perm[0]],
                                                    MatsubaraNumbers[where->FrequenciesPermutation.perm[1]],
                                                    MatsubaraNumbers[where->FrequenciesPermutation.perm[2]]
                                                   );
        return ((RealType)where->FrequenciesPermutation.sign)*Value;
    };
    return 0;
}

const std::vector<TwoParticleGFContainer::IndexCombination*>& TwoParticleGFContainer::getNonTrivialCombinations()
{
    return NonTrivialCombinations;
};

bool TwoParticleGFContainer::vanishes(const IndexCombination& in)
{
    return (mapNonTrivialCombinations.count(in) == 0);
}

long TwoParticleGFContainer::getNumberOfMatsubaras() const
{
    return NumberOfMatsubaras;
}

} // end of namespace Pomerol
