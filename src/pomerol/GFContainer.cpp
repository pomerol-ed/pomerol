#include "pomerol/GFContainer.h"

namespace Pomerol{

template<bool Complex>
GFContainer<Complex>::GFContainer (const IndexClassification<Complex>& IndexInfo,
                                   const StatesClassification<Complex>& S,
                                   const Hamiltonian<Complex> &H,
                                   const DensityMatrix<Complex> &DM,
                                   const FieldOperatorContainer<Complex>& Operators) :
    ContainerBase(this,IndexInfo),
    Thermal(DM), S(S), H(H), DM(DM), Operators(Operators)
{}

template<bool Complex>
void GFContainer<Complex>::prepareAll(const std::set<IndexCombination2>& InitialIndices)
{
    ContainerBase::fill(InitialIndices);
    for(auto iter = ContainerBase::ElementsMap.begin(); iter != ContainerBase::ElementsMap.end(); iter++)
        (iter->second)->prepare();
}

template<bool Complex>
void GFContainer<Complex>::computeAll()
{
    for(auto iter = ContainerBase::ElementsMap.begin(); iter != ContainerBase::ElementsMap.end(); iter++)
        (iter->second)->compute();
}

template<bool Complex>
auto GFContainer<Complex>::createElement(const IndexCombination2& Indices) const -> ElementT*
{
    return new ElementT(S,
                        H,
                        Operators.getAnnihilationOperator(Indices.Index1),
                        Operators.getCreationOperator(Indices.Index2),DM);
}

template class GFContainer<false>;
template class GFContainer<true>;

} // end of namespace Pomerol
