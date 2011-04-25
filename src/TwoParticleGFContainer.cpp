#include "TwoParticleGFContainer.h"

extern std::ostream& OUTPUT_STREAM;

const Permutation4 FourIndexContainer::TrivialOperatorPermutations[4] = { permutations4[0], permutations4[1], permutations4[6], permutations4[7] };

/*=========================================================================*/
TwoParticleGFContainer::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4)
{
    Indices[0]=cindex1;
    Indices[1]=cindex2;
    Indices[2]=cdagindex3;
    Indices[3]=cdagindex4;
}

bool TwoParticleGFContainer::IndexCombination::operator<(const TwoParticleGFContainer::IndexCombination& rhs) const
{
  return (Indices[0] < rhs.Indices[0]) || 
         (Indices[0] == rhs.Indices[0] && Indices[1] < rhs.Indices[1] ) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] < rhs.Indices[2]) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] < rhs.Indices[3]); 
}

bool TwoParticleGFContainer::IndexCombination::operator==(const TwoParticleGFContainer::IndexCombination& rhs) const
{
    return (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] == rhs.Indices[3]); 
}

bool TwoParticleGFContainer::IndexCombination::operator!=(const TwoParticleGFContainer::IndexCombination& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& output,const TwoParticleGFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

/*=========================================================================*/

TwoParticleGFContainer::FourIndexPermutation::FourIndexPermutation(const TwoParticleGFContainer::IndexCombination in, Permutation4 FrequenciesPermutation):
    to(TwoParticleGFContainer::IndexCombination(in.Indices[0],in.Indices[1],in.Indices[2],in.Indices[3])),FrequenciesPermutation(FrequenciesPermutation)
{}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,IndexClassification& IndexInfo, FieldOperatorContainer& Operators):
    S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators),NumberOfMatsubaras(0)
{
};


void TwoParticleGFContainer::readNonTrivialIndices(std::vector<IndexCombination*>& in)
{
NonTrivialCombinations=in;

for (std::vector<IndexCombination*>::const_iterator it1=in.begin(); it1!=in.end(); ++it1){
    DEBUG(**it1);
}
};

void TwoParticleGFContainer::prepare()
{
for (std::vector<IndexCombination*>::const_iterator it1=NonTrivialCombinations.begin(); it1!=NonTrivialCombinations.end(); ++it1){
    AnnihilationOperator &C1 = Operators.getAnnihilationOperator((*it1)->Indices[0]);
    AnnihilationOperator &C2 = Operators.getAnnihilationOperator((*it1)->Indices[1]);
    CreationOperator     &CX3 = Operators.getCreationOperator   ((*it1)->Indices[2]);
    CreationOperator     &CX4 = Operators.getCreationOperator   ((*it1)->Indices[3]);
    TwoParticleGF * temp2PGF = new TwoParticleGF(S,H,C1,C2,CX3,CX4,DM);
    temp2PGF->prepare();
    if (!temp2PGF->vanishes()){
        mapNonTrivialCombinations[**it1] = temp2PGF;
        INFO("TwoParticleGFContainer: assigned " << **it1 << " TwoParticleGF for calculation");
        addInnerPermutationsOfIndexCombination(*it1,*it1);
        for (std::list<IndexClassification::IndexPermutation>::const_iterator index_perm = IndexInfo.getEquivalentIndexPermutations().begin();
                index_perm!=IndexInfo.getEquivalentIndexPermutations().end(); ++index_perm){
                    IndexCombination* permuted = new IndexCombination(index_perm->Permutation[(*it1)->Indices[0]],
                                                                      index_perm->Permutation[(*it1)->Indices[1]],
                                                                      index_perm->Permutation[(*it1)->Indices[2]],
                                                                      index_perm->Permutation[(*it1)->Indices[3]]
                                                                     );
                    FourIndexPermutation* permuted_from = new FourIndexPermutation(**it1,TrivialOperatorPermutations[0]);
                    if (mapTrivialCombinations.count(*permuted)==0) 
                        mapTrivialCombinations[*permuted] = permuted_from;
                    else{
                        INFO_NONEWLINE("TwoParticleGFContainer: Detected frequency symmetry in " << *permuted << " : 1234=")
                        INFO_NONEWLINE((mapTrivialCombinations[*permuted]->FrequenciesPermutation.sign==-1?"-":""));
                        INFO_NONEWLINE(mapTrivialCombinations[*permuted]->FrequenciesPermutation.perm[0]+1);
                        INFO_NONEWLINE(mapTrivialCombinations[*permuted]->FrequenciesPermutation.perm[1]+1);
                        INFO_NONEWLINE(mapTrivialCombinations[*permuted]->FrequenciesPermutation.perm[2]+1);
                        INFO_NONEWLINE(mapTrivialCombinations[*permuted]->FrequenciesPermutation.perm[3]+1);
                        INFO("!");
                        };
                    TrivialCombinations.push_back(permuted);
                    INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << " = " << **it1 );
                    INFO(" with a frequencies permutation 1234->" << 1234 << " due to symmetry considerations");
                    addInnerPermutationsOfIndexCombination(permuted,*it1);
                };
        };
    };
};

void TwoParticleGFContainer::addInnerPermutationsOfIndexCombination(const IndexCombination *in, const IndexCombination *in_nontrivial)
{
    for (int perm = 1; perm < 4; ++perm){
        Permutation4 perm4 = TrivialOperatorPermutations[perm];
        IndexCombination *permuted = new IndexCombination(in->Indices[perm4.perm[0]],in->Indices[perm4.perm[1]],in->Indices[perm4.perm[2]],in->Indices[perm4.perm[3]]);
        if (*permuted!=*in){
            FourIndexPermutation* permuted_from = new FourIndexPermutation(*in_nontrivial,perm4);
            mapTrivialCombinations[*permuted] = permuted_from;
            TrivialCombinations.push_back(permuted);
            INFO_NONEWLINE("TwoParticleGFContainer: assigned " << *permuted << " = " << ((perm4.sign==-1)?"-":" ") << *in);
            INFO_NONEWLINE(" with a frequencies permutation 1234->" << perm4.perm[0]+1 << perm4.perm[1]+1 << perm4.perm[2]+1 << perm4.perm[3]+1);
            INFO(" due to internal indices permutation");
        };
    };
}

void TwoParticleGFContainer::compute(long NumberOfMatsubaras)
{
    this->NumberOfMatsubaras = NumberOfMatsubaras;
    for (std::map<IndexCombination,TwoParticleGF*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
        it1->second->compute(NumberOfMatsubaras);
        };
};


ComplexType TwoParticleGFContainer::operator()(const IndexCombination& in, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    if ( mapNonTrivialCombinations.count(in) > 0 ) return (*mapNonTrivialCombinations[in])(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
    if ( mapTrivialCombinations.count(in) > 0 ){
        FourIndexPermutation* where = mapTrivialCombinations[in];
        long MatsubaraNumbers[4] = {MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3,MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3};
        //DEBUG("Matsubaras: " << MatsubaraNumber1 <<"," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3);
       // DEBUG(in << " is equivalent to" << mapTrivialCombinations[in]->to);
        ComplexType Value = (*mapNonTrivialCombinations[where->to])(
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

const std::vector<TwoParticleGFContainer::IndexCombination*>& TwoParticleGFContainer::getTrivialCombinations()
{
    return TrivialCombinations;
};

bool TwoParticleGFContainer::vanishes(const IndexCombination& in)
{
    return ( mapNonTrivialCombinations.count(in) == 0  && mapTrivialCombinations.count(in) == 0);
}

long TwoParticleGFContainer::getNumberOfMatsubaras() const
{
    return NumberOfMatsubaras;
}

RealType TwoParticleGFContainer::getBeta() const
{
    return DM.getBeta();
}
