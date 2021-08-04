#include "pomerol/MonomialOperator.hpp"

#include <cstdlib>

namespace Pomerol {

void MonomialOperator::checkPrepared() const
{
    if(getStatus() < Prepared) {
        throw StatusMismatch("MonomialOperator is not prepared yet.");
    }
}

MonomialOperator::BlocksBimap const& MonomialOperator::getBlockMapping() const
{
    checkPrepared();
    return LeftRightBlocks;
}

void MonomialOperator::compute(const MPI_Comm& comm)
{
    checkPrepared();
    if(getStatus() >= Computed) return;

    std::size_t Size = parts.size();
    for(std::size_t BlockIn = 0; BlockIn < Size; BlockIn++){
        INFO_NONEWLINE( (int) ((1.0*BlockIn/Size) * 100 ) << "  " << std::flush);
        parts[BlockIn].compute();
    };

    setStatus(Computed);
}

MonomialOperatorPart& MonomialOperator::getPartFromRightIndex(BlockNumber out)
{
    checkPrepared();
    return parts[mapPartsFromRight.find(out)->second];
}

const MonomialOperatorPart& MonomialOperator::getPartFromRightIndex(BlockNumber out) const
{
    checkPrepared();
    return parts[mapPartsFromRight.find(out)->second];
}

MonomialOperatorPart& MonomialOperator::getPartFromLeftIndex(BlockNumber in)
{
    checkPrepared();
    return parts[mapPartsFromLeft.find(in)->second];
}

const MonomialOperatorPart& MonomialOperator::getPartFromLeftIndex(BlockNumber in) const
{
    checkPrepared();
    return parts[mapPartsFromLeft.find(in)->second];
}

BlockNumber MonomialOperator::getRightIndex(BlockNumber LeftIndex) const
{
    checkPrepared();
    BlocksBimap::left_const_iterator it =  LeftRightBlocks.left.find(LeftIndex);
    return it != LeftRightBlocks.left.end() ? it->second : INVALID_BLOCK_NUMBER;
}

BlockNumber MonomialOperator::getLeftIndex(BlockNumber RightIndex) const
{
    checkPrepared();
    BlocksBimap::right_const_iterator it =  LeftRightBlocks.right.find(RightIndex);
    return (it != LeftRightBlocks.right.end()) ? it->second : INVALID_BLOCK_NUMBER;
}

} // namespace Pomerol
