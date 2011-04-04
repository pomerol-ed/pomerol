#ifndef __CONTAINER2PGF__
#define __CONTAINER2PGF__

#include "config.h"
#include "TwoParticleGF.h"
#include <map>

//extern ostream OUTPUT_STREAM;

class TwoParticleGFContainer;


class TwoParticleGFContainer
{
public:
    struct IndexCombination;
    TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H);

    void readNonTrivialIndices(std::vector<IndexCombination*>&);
    void defineOperatorMaps();
  
    //void defineNonTrivialIndices();
private:
    StatesClassification S;
    Hamiltonian H;
    std::vector<IndexCombination*> NonTrivialCombinations;
    std::map<IndexCombination,TwoParticleGF*> NonTrivialValues;
    std::map<ParticleIndex,AnnihilationOperator*> mapAnnihilationOperators;
    std::map<ParticleIndex,CreationOperator*> mapCreationOperators;

};

struct TwoParticleGFContainer::IndexCombination
{
    ParticleIndex Indices[4];
    friend std::ostream& operator<<(std::ostream& output, const TwoParticleGFContainer::IndexCombination& out);
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};
#endif // endif :: #ifndef __CONTAINER2PGF__
