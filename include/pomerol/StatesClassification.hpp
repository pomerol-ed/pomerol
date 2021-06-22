/** \file include/pomerol/StatesClassification.h
** \brief Declaration of BlockNumber and StatesClassification classes
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_STATESCLASSIFICATION_H
#define POMEROL_INCLUDE_STATESCLASSIFICATION_H

#include "Misc.hpp"
#include "HilbertSpace.hpp"

#include <vector>

// TODO: Re-add functionality to classify many-body states by their quantum numbers

namespace Pomerol {

/** Index of a block (sector) within a many-body Hilbert space */
using BlockNumber = int;
/** All blocks with this number are treated as nonexistent */
constexpr BlockNumber INVALID_BLOCK_NUMBER = -1;

/** InnerQuantumState labels the states inside of the block of Fock States. Has no physical meaning. */
using InnerQuantumState = unsigned long;

/** This class handles all information about Fock states.
 *  It makes a classification of Fock states into blocks.
 */
class StatesClassification : public ComputableObject {

    /** A storage for all FockStates, each subvector correspond to the states
     *  which belong to a block with a given BlockNumber.
     */
    std::vector<std::vector<QuantumState>> StatesContainer;
    /** Index all states to belong to a block. */
    std::vector<BlockNumber> StateBlockIndex;

public:

    StatesClassification() = default;

    /** Perform a classification of all FockStates */
    // TODO: Optionally accept a list of integrals of motion
    // (Operators::expression<ScalarType, IndexTypes...> objects) and fill lists
    // of quantum numbers.
    template<typename ScalarType, typename... IndexTypes>
    void compute(HilbertSpace<ScalarType, IndexTypes...> const& HS) {
        if(Status == Computed) return;
        auto const& FullHilbertSpace = HS.getFullHilbertSpace();
        auto Dim = FullHilbertSpace.dim();
        if(HS.getStatus() == Computed) { // Multiple blocks revealed by HS
            InitMultipleBlocks(HS.getSpacePartition());
        } else { // Just one block
            InitSingleBlock(Dim);
        }
        Status = Computed;
    }

    /** Get total number of Quantum States */
    unsigned long getNumberOfStates() const { return StateBlockIndex.size(); }

    /** Returns total amount of non-vanishing blocks */
    BlockNumber getNumberOfBlocks() const { return StatesContainer.size(); }

    /**
     */
    unsigned long getBlockSize(BlockNumber in) const;

    /** get a vector of all FockStates with a given set of QuantumNumbers
     * \param[in] in A set of quantum numbers to get a vector of FockStates
     */
    const std::vector<QuantumState>& getFockStates(BlockNumber in) const;

    /** get a FockState, corresponding to an internal InnerQuantumState
     * \param[in] QuantumNumbers of block in which the InnerQuantumState is located
     * \param[in] m InnerQuantumState for which the correspondence is required
     */
    QuantumState getFockState(BlockNumber in, InnerQuantumState m) const;

    /** Returns BlockNumber of a given FockState
     * \param[in] in A FockState for which the BlockNumber is requested
     */
    BlockNumber getBlockNumber(QuantumState in) const;

    /** get InnerQuantumState of a given FockState. Since FockState is associated with
     * the Block number no explicit BlockNumber or QuantumNumbers is required
     * \param[in] state FockState for which the correspondence is required
     */
    InnerQuantumState getInnerState(QuantumState in) const;

private:

    void InitSingleBlock(unsigned long Dim);
    void InitMultipleBlocks(libcommute::space_partition const& partition);
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_STATESCLASSIFICATION_H
