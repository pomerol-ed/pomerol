#include "TwoParticleGFContainer.h"

extern output_handle OUT;
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
  return (Indices[0]*64+Indices[1]*16+Indices[2]*4+Indices[3] < rhs.Indices[0]*64+rhs.Indices[1]*16+rhs.Indices[2]*4+rhs.Indices[3]);
}


std::ostream& operator<<(std::ostream& output,const TwoParticleGFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM):S(S),H(H),DM(DM)
{
};


void TwoParticleGFContainer::readNonTrivialIndices(std::vector<IndexCombination*>& in)
{
NonTrivialCombinations=in;

for (std::vector<IndexCombination*>::const_iterator it1=in.begin(); it1!=in.end(); ++it1){
    DEBUG(**it1);
}
};

void TwoParticleGFContainer::defineFieldOperatorMaps()
{
for (std::vector<IndexCombination*>::const_iterator it1=NonTrivialCombinations.begin(); it1!=NonTrivialCombinations.end(); ++it1){
    for (unsigned short i=0;i<2;i++){
            ParticleIndex currentIndex=(*it1)->Indices[i];
            if (!mapAnnihilationOperators.count(currentIndex)){
                mapAnnihilationOperators[currentIndex]=new AnnihilationOperator(S,H,OUT,currentIndex);
                OUTPUT_STREAM << "Created c_" << currentIndex << endl; 
                }
            else OUTPUT_STREAM << "c_" << currentIndex << " exists." << endl;

        };
    for (unsigned short i=2;i<4;i++){
            ParticleIndex currentIndex=(*it1)->Indices[i];
            if (!mapCreationOperators.count(currentIndex)){
                OUTPUT_STREAM << "Created c^+_" << currentIndex << endl; 
                mapCreationOperators[currentIndex]=new CreationOperator(S,H,OUT,currentIndex);
                }
            else OUTPUT_STREAM << "c^+_" << currentIndex << " exists." << endl;

        };
 
    };
};

void TwoParticleGFContainer::computeFieldOperators()
{
for(std::map<ParticleIndex,AnnihilationOperator*>::iterator it1=mapAnnihilationOperators.begin();it1!=mapAnnihilationOperators.end();++it1){
    it1->second->prepare();
    it1->second->compute();
    };
for(std::map<ParticleIndex,CreationOperator*>::iterator it1=mapCreationOperators.begin();it1!=mapCreationOperators.end();++it1){
    it1->second->prepare();
    it1->second->compute();
    };
}

void TwoParticleGFContainer::prepareTwoParticleGFs()
{
for (std::vector<IndexCombination*>::const_iterator it1=NonTrivialCombinations.begin(); it1!=NonTrivialCombinations.end(); ++it1){
    AnnihilationOperator *C1 = mapAnnihilationOperators[(*it1)->Indices[0]];
    AnnihilationOperator *C2 = mapAnnihilationOperators[(*it1)->Indices[1]];
    CreationOperator     *CX3 = mapCreationOperators   [(*it1)->Indices[2]];
    CreationOperator     *CX4 = mapCreationOperators   [(*it1)->Indices[3]];
    mapNonTrivialCombinations[**it1]=new TwoParticleGF(S,H,*C1,*C2,*CX3,*CX4,DM);
    mapNonTrivialCombinations[**it1]->prepare();
    };
};


void TwoParticleGFContainer::computeTwoParticleGFs(long NumberOfMatsubaras)
{
for (std::map<IndexCombination,TwoParticleGF*>::iterator it1=mapNonTrivialCombinations.begin();it1!=mapNonTrivialCombinations.end();++it1){
    it1->second->compute(NumberOfMatsubaras);
    };
};


ComplexType TwoParticleGFContainer::operator()(IndexCombination& in, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
  return (*mapNonTrivialCombinations[in])(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}
