//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#include "StatesClassification.h"

namespace Pomerol{

bool BlockNumber::operator<(const BlockNumber& rhs) const {return number<rhs.number;}
bool BlockNumber::operator==(const BlockNumber& rhs) const {return number==rhs.number;}

//
// StatesClassification
//

StatesClassification::StatesClassification(const IndexClassification& IndexInfo, const Symmetrizer &Symm):IndexInfo(IndexInfo),Symm(Symm)
{
    IndexSize = IndexInfo.getIndexSize();
    StateSize = 1<<IndexSize;
}

void StatesClassification::compute()             
{
    std::vector<boost::shared_ptr<Operator> > sym_op = Symm.getOperations();
    int NOperations=sym_op.size();
    BlockNumber block_index=0;
    for (unsigned long FockStateIndex=0; FockStateIndex<StateSize; ++FockStateIndex) {
        FockState current_state(IndexSize,FockStateIndex);
        QuantumNumbers QNumbers(Symm.getQuantumNumbers());
        for (int n=0; n<NOperations; ++n) {
            RealType Value=sym_op[n]->getMatrixElement(current_state, current_state);
            QNumbers.set(n,Value);
        }
        std::map<QuantumNumbers, BlockNumber>::iterator map_pos=QuantumToBlock.find(QNumbers);
        if (map_pos==QuantumToBlock.end()) {
            DEBUG("Adding " << current_state << " to block " << block_index << " with QuantumNumbers " << QNumbers << ".");
            QuantumToBlock[QNumbers]=block_index;
            BlockToQuantum.insert(std::make_pair(block_index, QNumbers)); // Needed not to invoke an empty constructor.
            StatesContainer.push_back(std::vector<FockState>(0));
            StatesContainer[block_index].push_back(current_state);
            StateBlockIndex.push_back(block_index);
            block_index++;
            }
         else {
            DEBUG("Adding " << current_state << " to block " << map_pos->second << " with QuantumNumbers " << QNumbers << ".");
            StatesContainer[map_pos->second].push_back(current_state);
            StateBlockIndex.push_back(map_pos->second);
            };
        }
}

BlockNumber StatesClassification::getBlockNumber(QuantumNumbers in) const
{
    return (QuantumToBlock.count(in))?QuantumToBlock.find(in)->second:ERROR_BLOCK_NUMBER;
}

const unsigned long StatesClassification::getNumberOfStates() const
{
    return StateSize;
}

BlockNumber StatesClassification::getBlockNumber(FockState in) const
{
    if ( in.to_ulong() > StateSize ) { throw exWrongState(); };
    return StateBlockIndex[in.to_ulong()];
}

const InnerQuantumState StatesClassification::getInnerState(FockState state) const
{
  if ( state.to_ulong() > StateSize ) { throw (exWrongState()); return StateSize; };
  BlockNumber block = this->getBlockNumber(state);
  for (InnerQuantumState n=0; n<StatesContainer[block].size(); n++ )
    {    
        if ( StatesContainer[block][n]== state) { return n; } 
    }
 throw (exWrongState());
 return StateSize;
}

const std::vector<FockState>& StatesClassification::getFockStates( BlockNumber in ) const
{
    return StatesContainer[in];
}

const std::vector<FockState>& StatesClassification::getFockStates( QuantumNumbers in ) const
{
    std::map<QuantumNumbers,BlockNumber>::const_iterator it=QuantumToBlock.find(in);
    if (it != QuantumToBlock.end())
        return this->getFockStates(it->second);
    else
        throw (exWrongState());
}

const FockState StatesClassification::getFockState( BlockNumber in, InnerQuantumState m) const
{  
    if (int(in) < StatesContainer.size()) 
        if ( m < StatesContainer[in].size())
            return StatesContainer[in][m];
    ERROR("Couldn't find state numbered " << m << " in block " << in);
    throw (exWrongState());
    return ERROR_FOCK_STATE;
}

const FockState StatesClassification::getFockState( QuantumNumbers in, InnerQuantumState m) const
{
    return getFockState(getBlockNumber(in),m);
}


Symmetrizer::QuantumNumbers StatesClassification::getBlockInfo(BlockNumber in) const
{
    if (BlockToQuantum.count(in))
        return BlockToQuantum.find(in)->second;
    throw (exWrongState());
    return BlockToQuantum.find(0)->second;
}

BlockNumber StatesClassification::NumberOfBlocks() const
{
    return StatesContainer.size();
}

Symmetrizer::QuantumNumbers StatesClassification::getStateInfo(FockState in) const             
{
    return BlockToQuantum.find(getBlockNumber(in))->second;
}

const char* StatesClassification::exWrongState::what() const throw(){
    return "Wrong state";
};

} // end of namespace Pomerol

