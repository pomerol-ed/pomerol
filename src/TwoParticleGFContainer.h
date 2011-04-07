#ifndef __CONTAINER2PGF__
#define __CONTAINER2PGF__

#include "config.h"
#include "TwoParticleGF.h"
#include "FieldOperatorContainer.h"
#include <map>

class TwoParticleGFContainer;


class TwoParticleGFContainer
{
public:
    struct IndexCombination;
    TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM, BitClassification& IndexInfo, FieldOperatorContainer& Operators);

    void readNonTrivialIndices(std::vector<IndexCombination*>&);
    void prepare();
    void compute(long NumberOfMatsubaras);
    bool vanishes(IndexCombination&);
    long getNumberOfMatsubaras() const;

    ComplexType operator()(IndexCombination&, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3); 

    void dump();
  
    //void defineNonTrivialIndices();
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM; 
    BitClassification &IndexInfo;
    FieldOperatorContainer &Operators;
    std::vector<IndexCombination*> NonTrivialCombinations;
    std::map<IndexCombination,TwoParticleGF*> mapNonTrivialCombinations;
    long NumberOfMatsubaras;
};

struct TwoParticleGFContainer::IndexCombination
{
    ParticleIndex Indices[4];
    friend std::ostream& operator<<(std::ostream& output, const TwoParticleGFContainer::IndexCombination& out);
    bool operator<(const TwoParticleGFContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};
#endif // endif :: #ifndef __CONTAINER2PGF__
