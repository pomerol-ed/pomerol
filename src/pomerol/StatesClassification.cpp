#include "pomerol/StatesClassification.h"

namespace Pomerol{

bool BlockNumber::operator<(const BlockNumber& rhs) const {return number<rhs.number;}
bool BlockNumber::operator==(const BlockNumber& rhs) const {return number==rhs.number;}

//
// StatesClassification
//

StatesClassification::StatesClassification(const IndexClassification& IndexInfo, const Symmetrizer &Symm):
    ComputableObject(), 
    IndexInfo(IndexInfo),
    Symm(Symm)
{
}

void StatesClassification::compute()             
{
    if (Status>=Computed) return;
    IndexSize = IndexInfo.getIndexSize();
    StateSize = 1<<IndexSize;
    std::vector<boost::shared_ptr<Operator> > sym_op = Symm.getOperations();
    int NOperations=sym_op.size();
    BlockNumber block_index=0;
    for (QuantumState FockStateIndex=0; FockStateIndex<StateSize; ++FockStateIndex) {
        FockState current_state(IndexSize,FockStateIndex);
        QuantumNumbers QNumbers(Symm.getQuantumNumbers());
        for (int n=0; n<NOperations; ++n) {
            MelemType Value=sym_op[n]->getMatrixElement(current_state, current_state);
            QNumbers.set(n,Value);
        }
        std::map<QuantumNumbers, BlockNumber>::iterator map_pos=QuantumToBlock.find(QNumbers);
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

BlockNumber StatesClassification::getBlockNumber(QuantumNumbers in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return (QuantumToBlock.count(in))?QuantumToBlock.find(in)->second:ERROR_BLOCK_NUMBER;
}

const unsigned long StatesClassification::getNumberOfStates() const
{
    return StateSize;
}

BlockNumber StatesClassification::getBlockNumber(FockState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( in.to_ulong() > StateSize ) { throw exWrongState(); };
    return StateBlockIndex[in.to_ulong()];
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( in > StateSize ) { throw exWrongState(); };
    return StateBlockIndex[in];
}

const InnerQuantumState StatesClassification::getInnerState(FockState state) const
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

const InnerQuantumState StatesClassification::getInnerState(QuantumState state) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if ( state > StateSize ) { throw (exWrongState()); return StateSize; };
    return this->getInnerState(FockState(IndexSize, state));
}

const std::vector<FockState>& StatesClassification::getFockStates( BlockNumber in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return StatesContainer[in];
}

const std::vector<FockState>& StatesClassification::getFockStates( QuantumNumbers in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    std::map<QuantumNumbers,BlockNumber>::const_iterator it=QuantumToBlock.find(in);
    if (it != QuantumToBlock.end())
        return this->getFockStates(it->second);
    else
        throw (exWrongState());
}

const size_t StatesClassification::getBlockSize( BlockNumber in ) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return this->getFockStates(in).size();
}

const FockState StatesClassification::getFockState( BlockNumber in, InnerQuantumState m) const
{  
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if (int(in) < StatesContainer.size()) 
        if ( m < StatesContainer[in].size())
            return StatesContainer[in][m];
    ERROR("Couldn't find state numbered " << m << " in block " << in);
    throw (exWrongState());
    return ERROR_FOCK_STATE;
}

const FockState StatesClassification::getFockState( QuantumNumbers in, InnerQuantumState m) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    return getFockState(getBlockNumber(in),m);
}


Symmetrizer::QuantumNumbers StatesClassification::getQuantumNumbers(BlockNumber in) const
{
    if ( Status < Computed ) { ERROR("StatesClassification is not computed yet."); throw (exStatusMismatch()); };
    if (BlockToQuantum.count(in))
        return BlockToQuantum.find(in)->second;
    throw (exWrongState());
    return BlockToQuantum.find(0)->second;
}

BlockNumber StatesClassification::NumberOfBlocks() const
{
    return StatesContainer.size();
}

Symmetrizer::QuantumNumbers StatesClassification::getQuantumNumbers(FockState in) const             
{
    return BlockToQuantum.find(getBlockNumber(in))->second;
}

const char* StatesClassification::exWrongState::what() const throw(){
    return "Wrong state";
};

} // end of namespace Pomerol

