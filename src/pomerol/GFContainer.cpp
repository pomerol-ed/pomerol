#include "pomerol/GFContainer.hpp"

namespace Pomerol{

void GFContainer::prepareAll(const std::set<IndexCombination2>& InitialIndices)
{
    fill(InitialIndices);
    for(auto & el : ElementsMap)
        el.second->prepare();
}

void GFContainer::computeAll()
{
    for(auto & el : ElementsMap)
        el.second->compute();
}

std::shared_ptr<GreensFunction> GFContainer::createElement(const IndexCombination2& Indices) const
{
    return std::make_shared<GreensFunction>(
        S,
        H,
        Operators.getAnnihilationOperator(Indices.Index1),
        Operators.getCreationOperator(Indices.Index2),
        DM
    );
}

} // namespace Pomerol
