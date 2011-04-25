#ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
#define __INCLUDE_TWOPARTICLEGFCONTAINER_H

#include"Misc.h"
#include"ComputableObject.h"
#include"TwoParticleGF.h"
#include"FieldOperatorContainer.h"
class FourIndexContainer
{
protected:
    static const Permutation4 TrivialOperatorPermutations[];//[4] = { permutations4[0], permutations4[1], permutations4[6], permutations4[7] };
public:
    struct IndexCombination;
    struct FourIndexPermutation;
};

class TwoParticleGFContainer : public ComputableObject, public FourIndexContainer
{
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM; 
    IndexClassification &IndexInfo;
    FieldOperatorContainer &Operators;
    std::vector<IndexCombination*> NonTrivialCombinations;
    std::vector<IndexCombination*> TrivialCombinations;
    std::map<IndexCombination,TwoParticleGF*>         mapNonTrivialCombinations;
    std::map<IndexCombination,FourIndexPermutation*>  mapTrivialCombinations;
    long NumberOfMatsubaras;
public:
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
    const std::vector<TwoParticleGFContainer::IndexCombination*>& getTrivialCombinations();
private:
    void addInnerPermutationsOfIndexCombination(const IndexCombination *in, const IndexCombination *in_nontrivial);
};

struct TwoParticleGFContainer::IndexCombination
{
    ParticleIndex Indices[4];
    friend std::ostream& operator<<(std::ostream& output, const TwoParticleGFContainer::IndexCombination& out);
    bool operator< (const TwoParticleGFContainer::IndexCombination& rhs) const ;
    bool operator==(const TwoParticleGFContainer::IndexCombination& rhs) const ;
    bool operator!=(const TwoParticleGFContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};

struct TwoParticleGFContainer::FourIndexPermutation
{
    IndexCombination to;
    Permutation4 FrequenciesPermutation;
    FourIndexPermutation(const TwoParticleGFContainer::IndexCombination in, Permutation4 FrequenciesPermutation);
};
#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
