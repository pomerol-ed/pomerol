#include "pomerol/FieldOperatorContainer.h"

namespace Pomerol{

FieldOperatorContainer::FieldOperatorContainer(IndexClassification &IndexInfo, StatesClassification &S, const Hamiltonian &H, bool use_transpose) : 
    IndexInfo(IndexInfo), S(S), H(H), use_transpose(use_transpose)
{}

void FieldOperatorContainer::prepareAll(std::set<ParticleIndex> in)
{
    if (in.size() == 0) for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i) in.insert(i);
    for (std::set<ParticleIndex>::const_iterator it = in.begin(); it!=in.end(); it++)
        {
            ParticleIndex i = *it;
            CreationOperator *CX = new CreationOperator(IndexInfo, S,H,i);
            CX->prepare();
            mapCreationOperators[i] = CX;
            AnnihilationOperator *C = new AnnihilationOperator(IndexInfo, S,H,i);
            C->prepare();
            mapAnnihilationOperators[i] = C;
        }
}

void FieldOperatorContainer::computeAll()
{
    for (std::map <ParticleIndex, CreationOperator*>::iterator cdag_it = mapCreationOperators.begin(); cdag_it != mapCreationOperators.end(); ++cdag_it) {
        CreationOperator &cdag = *(cdag_it->second);
        cdag.compute();
        AnnihilationOperator &c = *mapAnnihilationOperators[cdag_it->first];

        FieldOperator::BlocksBimap cdag_block_map = cdag.getBlockMapping();
        // hack - copy transpose matrices into c
        for (FieldOperator::BlocksBimap::right_const_iterator cdag_map_it=cdag_block_map.right.begin(); cdag_map_it!=cdag_block_map.right.end(); cdag_map_it++) {
                c.getPartFromRightIndex(cdag_map_it->second).elementsRowMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getColMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).elementsColMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getRowMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).Status = ComputableObject::Computed;
                c.Status = ComputableObject::Computed;
            };
        };

// original
    //for (auto c : mapAnnihilationOperators) c.second->compute();
}

const CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in) const
{
    if (IndexInfo.checkIndex(in)){
        return *mapCreationOperators[in];
        }
    else
        throw (std::logic_error("No creation operator found."));
}

const AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in) const
{
    if (IndexInfo.checkIndex(in)){
        return *mapAnnihilationOperators[in];
        }
    else
        throw (std::logic_error("No annihilation operator found."));
}

} // end of namespace Pomerol
