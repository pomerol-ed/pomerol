#include "TwoParticleGFContainer.h"

extern std::ostream& OUTPUT_STREAM;

const Permutation4 FourIndexContainer::TrivialOperatorPermutations[4] = { permutations4[0], permutations4[1], permutations4[6], permutations4[7] };

/*=========================================================================*/
FourIndexContainer::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4)
{
    Indices[0]=cindex1;
    Indices[1]=cindex2;
    Indices[2]=cdagindex3;
    Indices[3]=cdagindex4;
}

bool FourIndexContainer::IndexCombination::operator<(const FourIndexContainer::IndexCombination& rhs) const
{
  return (Indices[0] < rhs.Indices[0]) || 
         (Indices[0] == rhs.Indices[0] && Indices[1] < rhs.Indices[1] ) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] < rhs.Indices[2]) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] < rhs.Indices[3]); 
}

bool FourIndexContainer::IndexCombination::operator==(const FourIndexContainer::IndexCombination& rhs) const
{
    return (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] == rhs.Indices[3]); 
}

bool FourIndexContainer::IndexCombination::operator!=(const FourIndexContainer::IndexCombination& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& output,const FourIndexContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

/*=========================================================================*/

TwoParticleGFContainer::Element::Element(TwoParticleGF* Computable2PGF, Permutation4 FrequenciesPermutation, bool isComputed):
    Computable2PGF(Computable2PGF),FrequenciesPermutation(FrequenciesPermutation),isComputed(isComputed)
{}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,IndexClassification& IndexInfo, FieldOperatorContainer& Operators):
    S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators),NumberOfMatsubaras(0)
{
};


void TwoParticleGFContainer::readInitialIndices(std::vector<IndexCombination*>& in)
{
InitialCombinations=in;

for (std::vector<IndexCombination*>::const_iterator it1=in.begin(); it1!=in.end(); ++it1){
}
};

void TwoParticleGFContainer::prepare()
{
if (Status < Prepared){
    for (std::vector<IndexCombination*>::const_iterator it1=InitialCombinations.begin(); it1!=InitialCombinations.end(); ++it1){
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
                    if (mapNonTrivialCombinations.count(*permuted)==0){
                        mapNonTrivialCombinations[*permuted] = permutedElement;
                        NonTrivialCombinations.push_back(permuted);
                        INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << " =  " << **it1 );
                        INFO(" with a frequencies permutation 1234->1234 due to symmetry considerations");
                        addInnerPermutationsOfIndexCombination(permuted);
                        }
                    else{
                        INFO_NONEWLINE("TwoParticleGFContainer: Detected frequency symmetry in " << *permuted << " : 1234=")
                        INFO_NONEWLINE((mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.sign==-1?"-":""));
                        INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[0]+1);
                        INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[1]+1);
                        INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[2]+1);
                        INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[3]+1);
                        INFO("!");
                        };
                };
            addInnerPermutationsOfIndexCombination(*it1);
                
        };
    };
    Status = Prepared;
    InitialCombinations.clear();
};
};

void TwoParticleGFContainer::addInnerPermutationsOfIndexCombination(const IndexCombination *in)
{
    IndexCombination *in_computed = new IndexCombination(
                                                         mapNonTrivialCombinations[*in]->Computable2PGF->getIndex(0), 
                                                         mapNonTrivialCombinations[*in]->Computable2PGF->getIndex(1),
                                                         mapNonTrivialCombinations[*in]->Computable2PGF->getIndex(2),
                                                         mapNonTrivialCombinations[*in]->Computable2PGF->getIndex(3)
                                                        );
    for (int perm = 1; perm < 4; ++perm){
        Permutation4 perm4 = TrivialOperatorPermutations[perm];
        IndexCombination *permuted = new IndexCombination(in->Indices[perm4.perm[0]],in->Indices[perm4.perm[1]],in->Indices[perm4.perm[2]],in->Indices[perm4.perm[3]]);
        if (*permuted!=*in){
            if (mapNonTrivialCombinations.count(*permuted)==0){
                Element* permutedElement = new Element(mapNonTrivialCombinations[*in_computed]->Computable2PGF,perm4,false);
                mapNonTrivialCombinations[*permuted] = permutedElement;
                NonTrivialCombinations.push_back(permuted);
                INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << " = " << ((perm4.sign==-1)?"-":" ") << *in_computed);
                INFO_NONEWLINE(" with a frequencies permutation 1234->" << perm4.perm[0]+1 << perm4.perm[1]+1 << perm4.perm[2]+1 << perm4.perm[3]+1);
                INFO(" due to internal indices permutation");
            }
            else if (mapNonTrivialCombinations[*permuted]->FrequenciesPermutation != TrivialOperatorPermutations[0]){
                // A thin moment. One checks that current permutation was not entered before i.e. by symmetry considerations or as an initial combination.
                // It might happen, that an obtained permutation is an exact copy not of IndexCombination *in, but rather the one
                // in NonTrivialCombinations, entered before. This 'if' operator is introduced to avoid stupid messages like 
                // "TwoParticleGFContainer: Detected frequency symmetry in ... : 1234=1234"
                INFO_NONEWLINE("TwoParticleGFContainer: Detected frequency symmetry in " << *permuted << " : 1234=");
                INFO_NONEWLINE((mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.sign==-1?"-":""));
                INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[0]+1);
                INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[1]+1);
                INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[2]+1);
                INFO_NONEWLINE(mapNonTrivialCombinations[*permuted]->FrequenciesPermutation.perm[3]+1);
                INFO("!");
                };

        };
    };
}

void TwoParticleGFContainer::compute(long NumberOfMatsubaras)
{
if (Status == Prepared){
    this->NumberOfMatsubaras = NumberOfMatsubaras;
    for (std::map<IndexCombination,Element*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
        it1->second->Computable2PGF->compute(NumberOfMatsubaras);
        };
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

RealType TwoParticleGFContainer::getBeta() const
{
    return DM.getBeta();
}
