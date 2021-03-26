#include "pomerol/StatesClassification.h"

#include <memory>

namespace Pomerol{

bool BlockNumber::operator<(const BlockNumber& rhs) const {return number<rhs.number;}
bool BlockNumber::operator==(const BlockNumber& rhs) const {return number==rhs.number;}

//
// StatesClassification
//

template<bool Complex>
StatesClassification<Complex>::StatesClassification(const IndexClassification<Complex>& IndexInfo, const Symmetrizer<Complex> &Symm):
    ComputableObject(),
    IndexInfo(IndexInfo),
    Symm(Symm)
{
}

template<bool Complex>
void StatesClassification<Complex>::compute()
{
    if (Status>=Computed) return;
    IndexSize = IndexInfo.getIndexSize();
    StateSize = 1<<IndexSize;
    std::vector<std::shared_ptr<Operator<Complex>> > sym_op = Symm.getOperations();
    int NOperations=sym_op.size();
    BlockNumber block_index=0;
    for (QuantumState FockStateIndex=0; FockStateIndex<StateSize; ++FockStateIndex) {
        FockState current_state(IndexSize,FockStateIndex);
        QuantumNumbers<Complex> QNumbers(Symm.getQuantumNumbers());
        for (int n=0; n<NOperations; ++n) {
            MelemType<Complex> Value=sym_op[n]->getMatrixElement(current_state, current_state);
            QNumbers.set(n,Value);
        }
        auto map_pos=QuantumToBlock.find(QNumbers);
        if (map_pos==QuantumToBlock.end()) {
//            DEBUG("Adding " << current_state << " to block " << block_index << " with QuantumNumbers " << QNumbers << ".");
            QuantumToBlock[QNumbers]=block_index;
            BlockToQuantum.insert(std::make_pair(block_index, QNumbers)); // Needed not to invoke an empty constructor.
            StatesContainer.push_back(std::vector<FockState>(0));
            StatesContainer[block_index].push_back(current_state);
            StateBlockIndex.push_back(block_index);
            block_index++;
            }
         else {
//            DEBUG("Adding " << current_state << " to block " << map_pos->second << " with QuantumNumbers " << QNumbers << ".");
            StatesContainer[map_pos->second].push_back(current_state);
            StateBlockIndex.push_back(map_pos->second);
            };
        }
    Status = Computed;
}

template<bool Complex>
BlockNumber StatesClassification<Complex>::getBlockNumber(QuantumNumbers<Complex> in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return (QuantumToBlock.count(in))?QuantumToBlock.find(in)->second:ERROR_BLOCK_NUMBER;
}

template<bool Complex>
const unsigned long StatesClassification<Complex>::getNumberOfStates() const
{
    return StateSize;
}

template<bool Complex>
BlockNumber StatesClassification<Complex>::getBlockNumber(FockState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( in.to_ulong() > StateSize ) { throw exWrongState(); };
    return StateBlockIndex[in.to_ulong()];
}

template<bool Complex>
BlockNumber StatesClassification<Complex>::getBlockNumber(QuantumState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( in > StateSize ) { throw exWrongState(); };
    return StateBlockIndex[in];
}

template<bool Complex>
const InnerQuantumState StatesClassification<Complex>::getInnerState(FockState state) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( state.to_ulong() > StateSize ) { throw (exWrongState()); return StateSize; };
    BlockNumber block = this->getBlockNumber(state);
    for (InnerQuantumState n=0; n<StatesContainer[block].size(); n++ )
      {
          if ( StatesContainer[block][n]== state) { return n; }
      }
    throw (exWrongState());
    return StateSize;
}

template<bool Complex>
const InnerQuantumState StatesClassification<Complex>::getInnerState(QuantumState state) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( state > StateSize ) { throw (exWrongState()); return StateSize; };
    return this->getInnerState(FockState(IndexSize, state));
}

template<bool Complex>
const std::vector<FockState>& StatesClassification<Complex>::getFockStates( BlockNumber in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return StatesContainer[in];
}

template<bool Complex>
const std::vector<FockState>& StatesClassification<Complex>::getFockStates( QuantumNumbers<Complex> in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    auto it = QuantumToBlock.find(in);
    if (it != QuantumToBlock.end())
        return this->getFockStates(it->second);
    else
        throw (exWrongState());
}

template<bool Complex>
const size_t StatesClassification<Complex>::getBlockSize( BlockNumber in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return this->getFockStates(in).size();
}

template<bool Complex>
const FockState StatesClassification<Complex>::getFockState( BlockNumber in, InnerQuantumState m) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if (int(in) < StatesContainer.size())
        if ( m < StatesContainer[in].size())
            return StatesContainer[in][m];
    ERROR("Couldn't find state numbered " << m << " in block " << in);
    throw (exWrongState());
    return ERROR_FOCK_STATE;
}

template<bool Complex>
const FockState StatesClassification<Complex>::getFockState( QuantumNumbers<Complex> in, InnerQuantumState m) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return getFockState(getBlockNumber(in),m);
}

template<bool Complex>
QuantumNumbers<Complex> StatesClassification<Complex>::getQuantumNumbers(BlockNumber in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if (BlockToQuantum.count(in))
        return BlockToQuantum.find(in)->second;
    throw (exWrongState());
    return BlockToQuantum.find(0)->second;
}

template<bool Complex>
BlockNumber StatesClassification<Complex>::NumberOfBlocks() const
{
    return StatesContainer.size();
}

template<bool Complex>
QuantumNumbers<Complex> StatesClassification<Complex>::getQuantumNumbers(FockState in) const
{
    return BlockToQuantum.find(getBlockNumber(in))->second;
}

template<bool Complex>
const char* StatesClassification<Complex>::exWrongState::what() const throw(){
    return "Wrong state";
};

template class StatesClassification<false>;
template class StatesClassification<true>;

} // end of namespace Pomerol

