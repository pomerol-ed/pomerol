#ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
#define __INCLUDE_TWOPARTICLEGFCONTAINER_H

#include"Misc.h"
#include"TwoParticleGF.h"
#include"FieldOperatorContainer.h"

class TwoParticleGFContainer
{
friend class Vertex4;

public:
    struct IndexCombination;
    TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM, IndexClassification& IndexInfo, FieldOperatorContainer& Operators);

    void readNonTrivialIndices(std::vector<IndexCombination*>&);
    void prepare();
    void compute(long NumberOfMatsubaras);
    bool vanishes(const IndexCombination&);
    long getNumberOfMatsubaras() const;

    ComplexType operator()(const IndexCombination&, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3); 

    void dump();
    RealType getBeta() const;
    const std::vector<TwoParticleGFContainer::IndexCombination*>& getNonTrivialCombinations();
  
    //void defineNonTrivialIndices();
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM; 
    IndexClassification &IndexInfo;
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
#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
