#include "pomerol/LibcommuteEigen.h"
#include "pomerol/MonomialOperator.h"

#include "mpi_dispatcher/mpi_skel.hpp"

namespace Pomerol {

MonomialOperator::BlocksBimap const& MonomialOperator::getBlockMapping() const
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }
    return LeftRightBlocks;
}

void MonomialOperator::compute(const MPI_Comm& comm)
{
    if (Status < Prepared) throw exStatusMismatch();
    if (Status >= Computed) return;

    size_t Size = parts.size();
    for (size_t BlockIn = 0; BlockIn < Size; BlockIn++){
        INFO_NONEWLINE( (int) ((1.0*BlockIn/Size) * 100 ) << "  " << std::flush);
        parts[BlockIn].compute();
    };
    Status = Computed;
}

MonomialOperatorPart& MonomialOperator::getPartFromRightIndex(BlockNumber out)
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }
    return parts[mapPartsFromRight.find(out)->second];
}

const MonomialOperatorPart& MonomialOperator::getPartFromRightIndex(BlockNumber out) const
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }
    return parts[mapPartsFromRight.find(out)->second];
}

MonomialOperatorPart& MonomialOperator::getPartFromLeftIndex(BlockNumber in)
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }
    return parts[mapPartsFromLeft.find(in)->second];
}

const MonomialOperatorPart& MonomialOperator::getPartFromLeftIndex(BlockNumber in) const
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }
    return parts[mapPartsFromLeft.find(in)->second];
}

BlockNumber MonomialOperator::getRightIndex(BlockNumber LeftIndex) const
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }

    BlocksBimap::left_const_iterator it =  LeftRightBlocks.left.find(LeftIndex);
    return (it != LeftRightBlocks.left.end()) ? it->second : INVALID_BLOCK_NUMBER;
}

BlockNumber MonomialOperator::getLeftIndex(BlockNumber RightIndex) const
{
    if (Status < Prepared) { ERROR("MonomialOperator is not prepared yet."); throw exStatusMismatch(); }

    BlocksBimap::right_const_iterator it =  LeftRightBlocks.right.find(RightIndex);
    return (it != LeftRightBlocks.right.end()) ? it->second : INVALID_BLOCK_NUMBER;
}

} // end of namespace Pomerol
