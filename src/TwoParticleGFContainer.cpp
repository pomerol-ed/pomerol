#include "TwoParticleGFContainer.h"

extern ostream& OUTPUT_STREAM;

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


std::ostream& operator<<(std::ostream& output,const TwoParticleGFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,BitClassification& IndexInfo, FieldOperatorContainer& Operators):
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
    if (!temp2PGF->vanishes()) mapNonTrivialCombinations[**it1] = temp2PGF;
    };
};


void TwoParticleGFContainer::compute(long NumberOfMatsubaras)
{
this->NumberOfMatsubaras = NumberOfMatsubaras;
for (std::map<IndexCombination,TwoParticleGF*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
    it1->second->compute(NumberOfMatsubaras);
    };
};


ComplexType TwoParticleGFContainer::operator()(const IndexCombination& in, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
if (!this->vanishes(in)) return (*mapNonTrivialCombinations[in])(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
else return 0;
}

const std::vector<TwoParticleGFContainer::IndexCombination*>& TwoParticleGFContainer::getNonTrivialCombinations()
{
    return NonTrivialCombinations;
};

bool TwoParticleGFContainer::vanishes(const IndexCombination& in)
{
return ( mapNonTrivialCombinations.count(in) == 0 );
}

long TwoParticleGFContainer::getNumberOfMatsubaras() const
{
    return NumberOfMatsubaras;
}

RealType TwoParticleGFContainer::getBeta() const
{
    return DM.getBeta();
}
