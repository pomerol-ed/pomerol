#ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
#define __INCLUDE_TWOPARTICLEGFCONTAINER_H

#include"Misc.h"
#include"ComputableObject.h"
#include"TwoParticleGF.h"
#include"FieldOperatorContainer.h"

class FourIndexContainer
{
protected:
    static const Permutation4 TrivialOperatorPermutations[];
public:
    struct IndexCombination;
};

class TwoParticleGFContainer : public ComputableObject, public FourIndexContainer
{
public:
    struct Element;
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM; 
    IndexClassification &IndexInfo;
    FieldOperatorContainer &Operators;
    std::vector<IndexCombination*> InitialCombinations;
    std::vector<IndexCombination*> NonTrivialCombinations;
    std::map<IndexCombination,Element*>  mapNonTrivialCombinations;
    long NumberOfMatsubaras;
    void defineInitialIndices();
public:
    TwoParticleGFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM, IndexClassification& IndexInfo, FieldOperatorContainer& Operators);

    void prepare();
    void compute(long NumberOfMatsubaras);
    bool vanishes(const IndexCombination&);
    long getNumberOfMatsubaras() const;

    ComplexType operator()(const IndexCombination&, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3); 

    void dump();
    RealType getBeta() const;
    const std::vector<FourIndexContainer::IndexCombination*>& getNonTrivialCombinations();
    const std::vector<FourIndexContainer::IndexCombination*>& getTrivialCombinations();
private:
    void addInnerPermutationsOfIndexCombination(const IndexCombination *in_nontrivial);
};

struct FourIndexContainer::IndexCombination
{
    ParticleIndex Indices[4];
    friend std::ostream& operator<<(std::ostream& output, const FourIndexContainer::IndexCombination& out);
    bool operator< (const FourIndexContainer::IndexCombination& rhs) const ;
    bool operator==(const FourIndexContainer::IndexCombination& rhs) const ;
    bool operator!=(const FourIndexContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};

struct TwoParticleGFContainer::Element
{
    TwoParticleGF* Computable2PGF;
    Permutation4 FrequenciesPermutation;
    bool isComputed;
    Element(TwoParticleGF* Computable2PGF, Permutation4 FrequenciesPermutation, bool isComputed);
    IndexCombination getLinkedIndices();
};
#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
