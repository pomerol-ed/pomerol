#ifndef __CONTAINER2PGF__
#define __CONTAINER2PGF__

#include "config.h"
#include "TwoParticleGF.h"
#include <map>

class TwoParticleGFContainer;


class TwoParticleGFContainer
{
public:
    struct IndexCombination;
    TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM);

    void readNonTrivialIndices(std::vector<IndexCombination*>&);
    void defineFieldOperatorMaps();
    void computeFieldOperators();
    void prepareTwoParticleGFs();
    void computeTwoParticleGFs(long NumberOfMatsubaras);

    ComplexType operator()(IndexCombination&, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3); 

    void dump();
  
    //void defineNonTrivialIndices();
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM; 
    std::vector<IndexCombination*> NonTrivialCombinations;
    std::map<IndexCombination,TwoParticleGF*> mapNonTrivialCombinations;
    std::map<ParticleIndex,AnnihilationOperator*> mapAnnihilationOperators;
    std::map<ParticleIndex,CreationOperator*> mapCreationOperators;

};

struct TwoParticleGFContainer::IndexCombination
{
    ParticleIndex Indices[4];
    friend std::ostream& operator<<(std::ostream& output, const TwoParticleGFContainer::IndexCombination& out);
    bool operator<(const TwoParticleGFContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};
#endif // endif :: #ifndef __CONTAINER2PGF__
