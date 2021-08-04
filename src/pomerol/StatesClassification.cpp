#include "pomerol/StatesClassification.hpp"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>

namespace Pomerol {

//
// StatesClassification
//

void StatesClassification::initSingleBlock(unsigned long Dim) {
    StateBlockIndex.resize(Dim, 0);
    StatesContainer.emplace_back(Dim, 0);
    std::iota(StatesContainer.back().begin(), StatesContainer.back().end(), 0);
}

void StatesClassification::initMultipleBlocks(libcommute::space_partition const& partition) {
    StateBlockIndex.resize(partition.dim());
    StatesContainer.resize(partition.n_subspaces(), std::vector<QuantumState>());
    foreach(partition, [this](QuantumState State, BlockNumber Block) {
        StateBlockIndex[State] = Block;
        StatesContainer[Block].push_back(State);
    });
}

void StatesClassification::checkComputed() const
{
    if(getStatus() < Computed) { throw StatusMismatch("StatesClassification is not computed yet."); }
}

unsigned long StatesClassification::getBlockSize(BlockNumber in) const
{
    checkComputed();
    return getFockStates(in).size();
}

const std::vector<QuantumState>& StatesClassification::getFockStates(BlockNumber in) const
{
    checkComputed();
    return StatesContainer[in];
}

QuantumState StatesClassification::getFockState(BlockNumber in, InnerQuantumState m) const
{
    checkComputed();
    if(int(in) < StatesContainer.size())
        if(m < StatesContainer[in].size())
            return StatesContainer[in][m];
    throw std::runtime_error("Wrong inner state " + std::to_string(m));
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const
{
    checkComputed();
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    return StateBlockIndex[in];
}

InnerQuantumState StatesClassification::getInnerState(QuantumState in) const
{
    checkComputed();
    if(in >= StateBlockIndex.size()) {
        throw std::runtime_error("Wrong state " + std::to_string(in));
    }
    BlockNumber n = StateBlockIndex[in];
    auto const& Block = StatesContainer[n];
    return std::distance(Block.begin(), std::find(Block.begin(), Block.end(), in));
}

} // namespace Pomerol
