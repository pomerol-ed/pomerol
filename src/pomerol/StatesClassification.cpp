#include "pomerol/StatesClassification.h"

#include <algorithm>
#include <string>
#include <numeric>
#include <stdexcept>

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
    if(getStatus() < Computed) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); }
    return getFockStates(in).size();
}

const std::vector<QuantumState>& StatesClassification::getFockStates(BlockNumber in) const
{
    if(getStatus() < Computed) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); }
    return StatesContainer[in];
}

QuantumState StatesClassification::getFockState(BlockNumber in, InnerQuantumState m) const
{
    if(getStatus() < Computed) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); }
    if(int(in) < StatesContainer.size())
        if(m < StatesContainer[in].size())
            return StatesContainer[in][m];
    ERROR("Couldn't find state numbered " << m << " in block " << in);
    throw std::runtime_error("Wrong inner state " + std::to_string(m));
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const
{
    if(getStatus() < Computed) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); }
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    return StateBlockIndex[in];
}

InnerQuantumState StatesClassification::getInnerState(QuantumState in) const
{
    if(getStatus() < Computed) { ERROR("StatesClassification is not computed yet."); throw exStatusMismatch(); }
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    BlockNumber n = StateBlockIndex[in];
    auto const& Block = StatesContainer[n];
    return std::distance(Block.begin(), std::find(Block.begin(), Block.end(), in));
}

} // namespace Pomerol
