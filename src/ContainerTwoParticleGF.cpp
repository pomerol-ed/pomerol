#include "ContainerTwoParticleGF.h"

extern output_handle OUT;
extern ostream& OUTPUT_STREAM;

TwoParticleGFContainer::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4)
{
    Indices[0]=cindex1;
    Indices[1]=cindex2;
    Indices[2]=cdagindex3;
    Indices[3]=cdagindex4;
}

std::ostream& operator<<(std::ostream& output,const TwoParticleGFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

/*=========================================================================*/

TwoParticleGFContainer::TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H):S(S),H(H)
{
};


void TwoParticleGFContainer::readNonTrivialIndices(std::vector<IndexCombination*>& in)
{
NonTrivialCombinations=in;

for (std::vector<IndexCombination*>::const_iterator it1=in.begin(); it1!=in.end(); ++it1){
    DEBUG(**it1);
}
};

void TwoParticleGFContainer::defineOperatorMaps()
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
