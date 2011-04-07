/** \file src/FieldOperatorContainer.cpp
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "FieldOperatorContainer.h"

FieldOperatorContainer::FieldOperatorContainer(StatesClassification &S, Hamiltonian &H, BitClassification &IndexInfo):S(S),H(H),IndexInfo(IndexInfo)
{
};

CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapCreationOperators.count(in)==0){
        CreationOperator *CX = new CreationOperator(S,H,in);
        DEBUG("Making Creation Operator_"<<in << std::endl);
        CX->prepare();
        CX->compute();
        mapCreationOperators[in] = CX;
        }
    else DEBUG("using Creation Operator_"<<in << std::endl);
    return *mapCreationOperators[in];
    }
else assert(0);
}

AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapAnnihilationOperators.count(in)==0){
        AnnihilationOperator *C = new AnnihilationOperator(S,H,in);
        DEBUG("Making Annihilation Operator_"<<in << std::endl);
        C->prepare();
        C->compute();
        mapAnnihilationOperators[in] = C;
        }
    else DEBUG("using Annihilation Operator_"<<in << std::endl);
    return *mapAnnihilationOperators[in];
    }
else assert(0);
}

