/** \file src/FieldOperatorContainer.cpp
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "FieldOperatorContainer.h"

FieldOperatorContainer::FieldOperatorContainer(StatesClassification &S, Hamiltonian &H, IndexClassification &IndexInfo):S(S),H(H),IndexInfo(IndexInfo)
{
};

CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapCreationOperators.count(in)==0){
        CreationOperator *CX = new CreationOperator(S,H,in);
        INFO("FieldOperatorContainer: Making Creation Operator_"<<in);
        CX->prepare();
        CX->compute();
        mapCreationOperators[in] = CX;
        };
    //else INFO("FieldOperatorContainer: Using already computed Creation Operator_"<< in);
    return *mapCreationOperators[in];
    }
else assert(0);
}

AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapAnnihilationOperators.count(in)==0){
        AnnihilationOperator *C = new AnnihilationOperator(S,H,in);
        INFO("FieldOperatorContainer: Making Annihilation Operator_"<<in);
        C->prepare();
        C->compute();
        mapAnnihilationOperators[in] = C;
        };
   // else INFO("FieldOperatorContainer: Using already computed Annihilation Operator_"<< in);
    return *mapAnnihilationOperators[in];
    }
else assert(0);
}

