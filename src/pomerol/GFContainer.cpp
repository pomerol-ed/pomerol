#include "pomerol/GFContainer.h"

namespace Pomerol{

GFContainer::GFContainer ( const IndexClassification& IndexInfo,
                           const StatesClassification& S,
                           const Hamiltonian &H, const DensityMatrix &DM,
                           const FieldOperatorContainer& Operators) :
    IndexContainer2<GreensFunction,GFContainer>(this,IndexInfo),
    Thermal(DM), S(S), H(H), DM(DM), Operators(Operators)
{}

void GFContainer::prepareAll(const std::set<IndexCombination2>& InitialIndices)
{
    fill(InitialIndices);
    for(std::map<IndexCombination2,GFPointer>::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        (iter->second)->prepare();
}

void GFContainer::computeAll()
{
    for(std::map<IndexCombination2,GFPointer>::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        (iter->second)->compute();
}

GreensFunction* GFContainer::createElement(const IndexCombination2& Indices) const
{
    return new GreensFunction(S,H, Operators.getAnnihilationOperator(Indices.Index1),
                                   Operators.getCreationOperator(Indices.Index2),DM);
}

} // end of namespace Pomerol
