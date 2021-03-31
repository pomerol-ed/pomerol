#include "pomerol/FieldOperator.h"

#include "mpi_dispatcher/mpi_skel.hpp"

namespace Pomerol{
FieldOperator::FieldOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    ComputableObject(), IndexInfo(IndexInfo), S(S), H(H), Index(Index)
{}

CreationOperator::CreationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,S,H,Index)
{
    O = new Pomerol::OperatorPresets::Cdag(Index);
}

AnnihilationOperator::AnnihilationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,S,H,Index)
{
    O = new Pomerol::OperatorPresets::C(Index);
}

FieldOperator::BlocksBimap const& FieldOperator::getBlockMapping() const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return LeftRightBlocks;
}

FieldOperatorPart& FieldOperator::getPartFromRightIndex(BlockNumber in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(in)->second];
}

FieldOperatorPart& FieldOperator::getPartFromRightIndex(const QuantumNumbers& in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(S.getBlockNumber(in))->second];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(BlockNumber in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(in)->second];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(const QuantumNumbers& in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(S.getBlockNumber(in))->second];
}

const std::vector<FieldOperatorPart*>& FieldOperator::getParts()
{
    return parts;
}

void FieldOperator::compute(const MPI_Comm& comm)
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

ParticleIndex FieldOperator::getIndex(void) const
{
    return Index;
}

void CreationOperator::prepare(void)
{
    if (Status >= Prepared) return;
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0; RightIndex<S.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        //DEBUG(RightIndex << "->" << LeftIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new CreationOperatorPart(IndexInfo, S,
                                    H.getPart(RightIndex),H.getPart(LeftIndex),Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("CreationOperator_" << Index <<": " << Size << " parts will be computed");
    Status = Prepared;
}

void AnnihilationOperator::prepare()
{
    if (Status >= Prepared) return;
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0;RightIndex<S.NumberOfBlocks();RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new AnnihilationOperatorPart(IndexInfo, S,
                                    H.getPart(RightIndex),H.getPart(LeftIndex), Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("AnnihilationOperator_" << Index <<": " << Size << " parts will be computed");
    Status = Prepared;
}

BlockNumber FieldOperator::getRightIndex(BlockNumber LeftIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }

    BlocksBimap::left_const_iterator it =  LeftRightBlocks.left.find(LeftIndex);
    return (it != LeftRightBlocks.left.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

BlockNumber FieldOperator::getLeftIndex(BlockNumber RightIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }

    BlocksBimap::right_const_iterator it =  LeftRightBlocks.right.find(RightIndex);
    return (it != LeftRightBlocks.right.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

BlockNumber FieldOperator::mapsTo(BlockNumber RightIndex) const
{
    bool found=false;
    std::map<FockState, ComplexType> result;
    const std::vector<FockState> &states=S.getFockStates(RightIndex);
    for (std::vector<FockState>::const_iterator state_it=states.begin(); state_it!=states.end() && !found; state_it++) {
        result = O->actRight(*state_it);
        found = (result.size()>0);
        }
    return (found)?S.getBlockNumber(result.begin()->first):ERROR_BLOCK_NUMBER;
}

QuantumNumbers FieldOperator::mapsTo(const QuantumNumbers& in) const
{
    BlockNumber out = this->mapsTo(S.getBlockNumber(in));
    if ( out == ERROR_BLOCK_NUMBER) throw (QuantumNumbers::exWrongNumbers());
    return S.getQuantumNumbers(out);
}

QuadraticOperator::QuadraticOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index1, ParticleIndex Index2) :
        FieldOperator(IndexInfo,S,H,9999), Index1(Index1), Index2(Index2)
        // Index=9999 dummy
{
    O = new Pomerol::OperatorPresets::N_offdiag(Index1, Index2);
}

void QuadraticOperator::prepare(void)
{
    if (Status >= Prepared) return;
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0; RightIndex<S.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new QuadraticOperatorPart(IndexInfo, S,
                    H.getPart(RightIndex), H.getPart(LeftIndex), Index1, Index2);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("QuadraticOperator_" << Index1 << "_" << Index2 <<": " << Size << " parts will be computed");
    Status = Prepared;
}

} // end of namespace Pomerol
