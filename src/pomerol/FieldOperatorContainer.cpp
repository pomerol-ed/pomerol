#include "pomerol/FieldOperatorContainer.h"

namespace Pomerol{

template<bool Complex>
FieldOperatorContainer<Complex>::FieldOperatorContainer(IndexClassification<Complex> &IndexInfo,
                                                        StatesClassification<Complex> &S,
                                                        const Hamiltonian<Complex> &H,
                                                        bool use_transpose) :
    IndexInfo(IndexInfo), S(S), H(H), use_transpose(use_transpose)
{}

template<bool Complex>
void FieldOperatorContainer<Complex>::prepareAll(std::set<ParticleIndex> in)
{
    if (in.size() == 0) for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i) in.insert(i);
    for (std::set<ParticleIndex>::const_iterator it = in.begin(); it!=in.end(); it++)
        {
            ParticleIndex i = *it;
            CreationOperator<Complex> *CX = new CreationOperator<Complex>(IndexInfo, S,H,i);
            CX->prepare();
            mapCreationOperators[i] = CX;
            AnnihilationOperator<Complex> *C = new AnnihilationOperator<Complex>(IndexInfo, S,H,i);
            C->prepare();
            mapAnnihilationOperators[i] = C;
        }
}

template<bool Complex>
void FieldOperatorContainer<Complex>::computeAll()
{
    for (auto cdag_it = mapCreationOperators.begin(); cdag_it != mapCreationOperators.end(); ++cdag_it) {
        CreationOperator<Complex> &cdag = *(cdag_it->second);
        cdag.compute();
        AnnihilationOperator<Complex> &c = *mapAnnihilationOperators[cdag_it->first];

        typename FieldOperator<Complex>::BlocksBimap cdag_block_map = cdag.getBlockMapping();
        // hack - copy transpose matrices into c
        for (auto cdag_map_it=cdag_block_map.right.begin(); cdag_map_it!=cdag_block_map.right.end(); cdag_map_it++) {
                c.getPartFromRightIndex(cdag_map_it->second).elementsRowMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getColMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).elementsColMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getRowMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).Status = ComputableObject::Computed;
                c.Status = ComputableObject::Computed;
            };
        };
}

template<bool Complex>
auto FieldOperatorContainer<Complex>::getCreationOperator(ParticleIndex in) const -> const CreationOperator<Complex>&
{
    if (IndexInfo.checkIndex(in)){
        return *mapCreationOperators[in];
        }
    else
        throw (std::logic_error("No creation operator found."));
}

template<bool Complex>
auto FieldOperatorContainer<Complex>::getAnnihilationOperator(ParticleIndex in) const -> const AnnihilationOperator<Complex>&
{
    if (IndexInfo.checkIndex(in)){
        return *mapAnnihilationOperators[in];
        }
    else
        throw (std::logic_error("No annihilation operator found."));
}

template class FieldOperatorContainer<false>;
template class FieldOperatorContainer<true>;

} // end of namespace Pomerol
