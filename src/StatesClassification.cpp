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
    std::list<boost::shared_ptr<Operator> > symmetric_operations_list = Symm.getOperations();
    DEBUG(symmetric_operations_list.size());
    for (unsigned long FockStateIndex=0; FockStateIndex<StateSize; ++FockStateIndex) {
        FockState current_state(IndexSize,FockStateIndex);
        DEBUG(current_state);
        for (std::list<boost::shared_ptr<Operator> >::const_iterator operations_it=symmetric_operations_list.begin(); operations_it!=symmetric_operations_list.end(); operations_it++) { 
            DEBUG((*operations_it)->getMatrixElement(current_state, current_state));   
        }
    }
}

/*
const InnerFockState StatesClassification::getInnerState(FockState state) const
{
  if ( state.to_ulong() > StateSize ) { throw exWrongState(); return StateSize; };
  for (InnerFockState n=0; n<(*this).getFockStates((*this).getStateInfo(state)).size(); n++ )
    {    
        if ( (*this).getFockState((*this).getStateInfo(state),n)== state) { return n; }   //get ST
    }
 return StateSize;
}

BlockNumber StatesClassification::getBlockNumber(QuantumNumbers in) const
{
    return (QuantumToBlock.count(in))?QuantumToBlock.find(in)->second:ERROR_BLOCK_NUMBER;
}

QuantumNumbers StatesClassification::getBlockInfo(BlockNumber in) const
{
      return (BlockToQuantum.count(in))?BlockToQuantum.find(in)->second:ERROR_QUANTUM_NUMBERS;
}

BlockNumber StatesClassification::NumberOfBlocks() const
{
}
const FockState StatesClassification::getNumberOfStates() const
{
    return StateSize;
}


QuantumNumbers StatesClassification::getStateInfo(FockState in) const             
{
}

BlockNumber StatesClassification::getBlockNumber(FockState in) const
{
//    return (*this).getBlockNumber((*this).getStateInfo(in));
}

const char* StatesClassification::exWrongState::what() const throw(){
    return "Wrong state";
};
*/

} // end of namespace Pomerol

