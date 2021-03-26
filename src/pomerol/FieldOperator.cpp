#include "pomerol/FieldOperator.h"

#include "mpi_dispatcher/mpi_skel.hpp"

namespace Pomerol{

template<bool Complex>
FieldOperator<Complex>::FieldOperator(const IndexClassification<Complex> &IndexInfo,
                                      const StatesClassification<Complex> &S,
                                      const Hamiltonian<Complex> &H,
                                      ParticleIndex Index) :
    ComputableObject(), IndexInfo(IndexInfo), S(S), H(H), Index(Index)
{}

template<bool Complex>
CreationOperator<Complex>::CreationOperator(const IndexClassification<Complex> &IndexInfo,
                                            const StatesClassification<Complex> &S,
                                            const Hamiltonian<Complex> &H, ParticleIndex Index) :
    FieldOperator<Complex>(IndexInfo,S,H,Index)
{
    Base::O = new Pomerol::OperatorPresets::Cdag<Complex>(Index);
}

template<bool Complex>
AnnihilationOperator<Complex>::AnnihilationOperator(const IndexClassification<Complex> &IndexInfo,
                                                    const StatesClassification<Complex> &S,
                                                    const Hamiltonian<Complex> &H, ParticleIndex Index) :
    FieldOperator<Complex>(IndexInfo,S,H,Index)
{
    Base::O = new Pomerol::OperatorPresets::C<Complex>(Index);
}

template<bool Complex>
auto FieldOperator<Complex>::getBlockMapping() const -> FieldOperator<Complex>::BlocksBimap const&
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return LeftRightBlocks;
}

template<bool Complex>
auto FieldOperator<Complex>::getPartFromRightIndex(BlockNumber in) const -> PartT&
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(in)->second];
}

template<bool Complex>
auto FieldOperator<Complex>::getPartFromRightIndex(const QuantumNumbers<Complex>& in) const -> PartT&
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(S.getBlockNumber(in))->second];
}

template<bool Complex>
auto FieldOperator<Complex>::getPartFromLeftIndex(BlockNumber in) const -> PartT&
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(in)->second];
}

template<bool Complex>
auto FieldOperator<Complex>::getPartFromLeftIndex(const QuantumNumbers<Complex>& in) const -> PartT&
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(S.getBlockNumber(in))->second];
}

template<bool Complex>
auto FieldOperator<Complex>::getParts() -> const std::vector<PartT*>&
{
    return parts;
}

template<bool Complex>
void FieldOperator<Complex>::compute(const MPI_Comm& comm)
{
    if (Status < Prepared) throw (exStatusMismatch());
    if (Status >= Computed) return;

    if (!pMPI::rank(comm))
      INFO_NONEWLINE("Computing " << *O << " in eigenbasis of the Hamiltonian: ");

    size_t Size = parts.size();
    for (size_t BlockIn = 0; BlockIn < Size; BlockIn++){
        INFO_NONEWLINE( (int) ((1.0*BlockIn/Size) * 100 ) << "  " << std::flush);
        parts[BlockIn]->compute();
    };
    INFO("");
    Status = Computed;
}

template<bool Complex>
ParticleIndex FieldOperator<Complex>::getIndex(void) const
{
    return Index;
}

template<bool Complex>
void CreationOperator<Complex>::prepare(void)
{
    if (Base::Status >= ComputableObject::Prepared) return;
    size_t Size = Base::parts.size();
    for (BlockNumber RightIndex=0; RightIndex<Base::S.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = Base::mapsTo(RightIndex);
        //DEBUG(RightIndex << "->" << LeftIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart<Complex> *Part = new CreationOperatorPart<Complex>(
                Base::IndexInfo,
                Base::S,
                Base::H.getPart(RightIndex),
                Base::H.getPart(LeftIndex),
                Base::Index
            );
            Base::parts.push_back(Part);
            Base::mapPartsFromRight[RightIndex]=Size;
            Base::mapPartsFromLeft[LeftIndex]=Size;
            Base::LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("CreationOperator_" << Base::Index <<": " << Size << " parts will be computed");
    Base::Status = ComputableObject::Prepared;
}

template<bool Complex>
void AnnihilationOperator<Complex>::prepare()
{
    if (Base::Status >= ComputableObject::Prepared) return;
    size_t Size = Base::parts.size();
    for (BlockNumber RightIndex=0;RightIndex<Base::S.NumberOfBlocks();RightIndex++){
        BlockNumber LeftIndex = Base::mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart<Complex> *Part = new AnnihilationOperatorPart<Complex>(
                Base::IndexInfo,
                Base::S,
                Base::H.getPart(RightIndex),
                Base::H.getPart(LeftIndex),
                Base::Index
            );
            Base::parts.push_back(Part);
            Base::mapPartsFromRight[RightIndex]=Size;
            Base::mapPartsFromLeft[LeftIndex]=Size;
            Base::LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("AnnihilationOperator_" << Base::Index <<": " << Size << " parts will be computed");
    Base::Status = ComputableObject::Prepared;
}

template<bool Complex>
BlockNumber FieldOperator<Complex>::getRightIndex(BlockNumber LeftIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }

    BlocksBimap::left_const_iterator it =  LeftRightBlocks.left.find(LeftIndex);
    return (it != LeftRightBlocks.left.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

template<bool Complex>
BlockNumber FieldOperator<Complex>::getLeftIndex(BlockNumber RightIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }

    BlocksBimap::right_const_iterator it =  LeftRightBlocks.right.find(RightIndex);
    return (it != LeftRightBlocks.right.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

template<bool Complex>
BlockNumber FieldOperator<Complex>::mapsTo(BlockNumber RightIndex) const
{
    bool found=false;
    std::map<FockState, MelemType<Complex>> result;
    const std::vector<FockState> &states=S.getFockStates(RightIndex);
    for (std::vector<FockState>::const_iterator state_it=states.begin(); state_it!=states.end() && !found; state_it++) {
        result = O->actRight(*state_it);
        found = (result.size()>0);
        }
    return (found)?S.getBlockNumber(result.begin()->first):ERROR_BLOCK_NUMBER;
}

template<bool Complex>
QuantumNumbers<Complex> FieldOperator<Complex>::mapsTo(const QuantumNumbers<Complex>& in) const
{
    BlockNumber out = this->mapsTo(S.getBlockNumber(in));
    if ( out == ERROR_BLOCK_NUMBER) throw (typename QuantumNumbers<Complex>::exWrongNumbers());
    return S.getQuantumNumbers(out);
}

template<bool Complex>
QuadraticOperator<Complex>::QuadraticOperator(const IndexClassification<Complex> &IndexInfo,
                                              const StatesClassification<Complex> &S,
                                              const Hamiltonian<Complex> &H,
                                              ParticleIndex Index1, ParticleIndex Index2) :
        FieldOperator<Complex>(IndexInfo,S,H,9999), Index1(Index1), Index2(Index2)
        // Index=9999 dummy
{
    Base::O = new Pomerol::OperatorPresets::N_offdiag<Complex>(Index1, Index2);
}

template<bool Complex>
void QuadraticOperator<Complex>::prepare(void)
{
    if (Base::Status >= ComputableObject::Prepared) return;
    size_t Size = Base::parts.size();
    for (BlockNumber RightIndex=0; RightIndex<Base::S.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = Base::mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart<Complex> *Part = new QuadraticOperatorPart<Complex>(
                Base::IndexInfo,
                Base::S,
                Base::H.getPart(RightIndex),
                Base::H.getPart(LeftIndex),
                Index1,
                Index2
            );
            Base::parts.push_back(Part);
            Base::mapPartsFromRight[RightIndex]=Size;
            Base::mapPartsFromLeft[LeftIndex]=Size;
            Base::LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("QuadraticOperator_" << Index1 << "_" << Index2 <<": " << Size << " parts will be computed");
    Base::Status = ComputableObject::Prepared;
}

// Explicit instantiations: Real case

template class FieldOperator<false>;
template class AnnihilationOperator<false>;
template class CreationOperator<false>;
template class QuadraticOperator<false>;

// Explicit instantiations: Complex case

template class FieldOperator<true>;
template class AnnihilationOperator<true>;
template class CreationOperator<true>;
template class QuadraticOperator<true>;

} // end of namespace Pomerol
