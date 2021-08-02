#include "pomerol/FieldOperatorContainer.hpp"

#include <stdexcept>

namespace Pomerol {

void FieldOperatorContainer::computeAll()
{
    for(auto & cdag_p : mapCreationOperators) {
        auto & cdag = cdag_p.second;
        cdag.compute();

        auto & c = mapAnnihilationOperators.find(cdag_p.first)->second;

        const auto & cdag_block_map = cdag_p.second.getBlockMapping();
        for(auto cdag_map_it = cdag_block_map.right.begin(); cdag_map_it != cdag_block_map.right.end(); ++cdag_map_it) {
            auto & cPart = c.getPartFromRightIndex(cdag_map_it->second);
            auto & cdagPart = cdag.getPartFromRightIndex(cdag_map_it->first);
            cPart.setFromAdjoint(cdagPart);
        }
        c.setStatus(ComputableObject::Computed);
    }
}

const CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in) const
{
    auto it = mapCreationOperators.find(in);
    if(it == mapCreationOperators.end())
        throw std::logic_error("No creation operator found.");
    else
        return it->second;
}

const AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in) const
{
    auto it = mapAnnihilationOperators.find(in);
    if(it == mapAnnihilationOperators.end())
        throw std::logic_error("No annihilation operator found.");
    else
        return it->second;
}

} // namespace Pomerol
