#include "pomerol/StatesClassification.h"

#include <algorithm>
#include <numeric>

namespace Pomerol{

//
// StatesClassification
//

void StatesClassification::InitSingleBlock(unsigned long Dim) {
    StateBlockIndex.resize(Dim, 0);
    StatesContainer.emplace_back(Dim, 0);
    std::iota(StatesContainer.back().begin(), StatesContainer.back().end(), 0);
}

void StatesClassification::InitMultipleBlocks(libcommute::space_partition const& partition) {
    StateBlockIndex.resize(partition.dim());
    StatesContainer.resize(partition.n_subspaces(), std::vector<QuantumState>());
    foreach(partition, [this](QuantumState State, BlockNumber Block) {
        StateBlockIndex[State] = Block;
        StatesContainer[Block].push_back(State);
    });
}

unsigned long StatesClassification::getBlockSize(BlockNumber in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return getFockStates(in).size();
}

const std::vector<QuantumState>& StatesClassification::getFockStates(BlockNumber in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return StatesContainer[in];
}

QuantumState StatesClassification::getFockState(BlockNumber in, InnerQuantumState m) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if (int(in) < StatesContainer.size())
        if ( m < StatesContainer[in].size())
            return StatesContainer[in][m];
    ERROR("Couldn't find state numbered " << m << " in block " << in);
    throw exWrongState(m); // FIXME
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); };
    if ( in >= StateBlockIndex.size() ) { throw exWrongState(in); };
    return StateBlockIndex[in];
}

InnerQuantumState StatesClassification::getInnerState(QuantumState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( in >= StateBlockIndex.size() ) { throw exWrongState(in); }
    BlockNumber n = StateBlockIndex[in];
    auto const& Block = StatesContainer[n];
    return std::distance(Block.begin(), std::find(Block.begin(), Block.end(), in));
}

} // end of namespace Pomerol
