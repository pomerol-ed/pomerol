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


const InnerQuantumState StatesClassification::getInnerState(QuantumState state) const
{
  InnerQuantumState ST=N_state;                // "state" in part of Hamilt            
  for (unsigned int n=0; n<(*this).getQuantumStates((*this).getStateInfo(state)).size(); n++ )
    {    
        if ( (*this).getQuantumState((*this).getStateInfo(state),n)== state) ST=n;    //get ST
    }
 return ST;
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
const QuantumState StatesClassification::getNumberOfStates() const
{
    return StateSize;
}

void StatesClassification::compute()             
{
}

QuantumNumbers StatesClassification::getStateInfo(QuantumState in) const             
{
}

BlockNumber StatesClassification::getBlockNumber(QuantumState in) const
{
//    return (*this).getBlockNumber((*this).getStateInfo(in));
}

} // end of namespace Pomerol

