/** \file src/GFContainer.h
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Andrey Antipov (antipov@ct-qmc.org)
*/


#ifndef ___GF_CONTAINER___
#define ___GF_CONTAINER___

#include "config.h"
#include "GreensFunction.h"
#include "FieldOperatorContainer.h"
#include <map>


class GFContainer
{

public:
    struct IndexCombination; 
    GFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM, BitClassification& IndexInfo, FieldOperatorContainer& Operators);
    void prepare();
    void compute();
    MatrixType& operator()(long MatsubaraNumber);
    ComplexType operator()(ParticleIndex i, ParticleIndex j, long MatsubaraNumber);
    void dumpToPlainText(long wn);
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM;
    BitClassification &IndexInfo;
    FieldOperatorContainer &Operators;

    std::map<IndexCombination, GreensFunction*> mapGreensFunctions;
};

struct GFContainer::IndexCombination
{
    ParticleIndex Indices[2];
    friend std::ostream& operator<<(std::ostream& output, const GFContainer::IndexCombination& out);
    bool operator<(const GFContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cdagindex2);
};

#endif // endif :: #ifndef ___GF_CONTAINER___
