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


#include "FieldOperator.h"
extern std::ostream& OUTPUT_STREAM;

FieldOperator::FieldOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index) : ComputableObject(),
    IndexInfo(IndexInfo), System(System), H(H), Index(Index)
{}

CreationOperator::CreationOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index) : 
    FieldOperator(IndexInfo,System,H,Index)
{}

AnnihilationOperator::AnnihilationOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index) : 
    FieldOperator(IndexInfo,System,H,Index)
{}

std::list<BlockMapping>& FieldOperator::getNonTrivialIndices()
{
    return LeftRightIndices;
};

FieldOperatorPart& FieldOperator::getPartFromRightIndex(BlockNumber in)
{
  return *Data[mapPartsFromRight[in]];
}

FieldOperatorPart& FieldOperator::getPartFromRightIndex(QuantumNumbers in)
{
  return *Data[mapPartsFromRight[System.getBlockNumber(in)]];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(BlockNumber in)
{
  return *Data[mapPartsFromLeft[in]];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(QuantumNumbers in)
{
  return *Data[mapPartsFromLeft[System.getBlockNumber(in)]];
}

void FieldOperator::compute()
{
if (Status < Computed ){
    size_t size = Data.size();
    INFO_NONEWLINE("FieldOperator_" << Index << ", computing: ")
    for (unsigned int b_in=0;b_in<size;b_in++){
        INFO_NONEWLINE( (int) ((1.0*b_in/size) * 100 ) << "  " << std::flush);
        Data[b_in]->compute();
        };
    INFO("");
    Status=Computed;
    };
}

ParticleIndex FieldOperator::getIndex() const
{ 
    return Index;
}

void CreationOperator::prepare()
{
if (Status < Prepared){
    size_t size = Data.size();
    for (BlockNumber RightIndex=0;RightIndex<System.NumberOfBlocks();RightIndex++){
            BlockNumber LeftIndex = this->mapsTo(RightIndex);
            if (LeftIndex.isCorrect()){
                FieldOperatorPart *Part = new CreationOperatorPart(IndexInfo, System,H.part(RightIndex),H.part(LeftIndex), Index);
                Data.push_back(Part);
                //OUTPUT_STREAM << "Entering creation operator part " << System.getBlockInfo(RightIndex) << "->" << System.getBlockInfo(LeftIndex) <<std::endl; 
                mapPartsFromRight[RightIndex]=size;
                mapPartsFromLeft[LeftIndex]=size;
                LeftRightIndices.push_back(BlockMapping(LeftIndex,RightIndex));
	            mapRightToLeftIndex[RightIndex]=LeftIndex;
	            mapLeftToRightIndex[LeftIndex]=RightIndex;
                size++;
                }    
            }
    Status=Prepared;
    INFO("CreationOperator_" << Index <<": " << size << " parts will be computed");
    };
}

void AnnihilationOperator::prepare()
{
if (Status < Prepared){
    size_t size = Data.size();
    for (BlockNumber RightIndex=0;RightIndex<System.NumberOfBlocks();RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new AnnihilationOperatorPart(IndexInfo, System,H.part(RightIndex),H.part(LeftIndex), Index);
            Data.push_back(Part);
            //OUTPUT_STREAM << "Entering annihilation operator part " << System.getBlockInfo(RightIndex) << "->" << System.getBlockInfo(LeftIndex) << std::endl; 
            mapPartsFromRight[RightIndex]=size;
            mapPartsFromLeft[LeftIndex]=size;
	        mapRightToLeftIndex[RightIndex]=LeftIndex;
	        mapLeftToRightIndex[LeftIndex]=RightIndex;
            LeftRightIndices.push_back(BlockMapping(LeftIndex,RightIndex));
            size++;
            }    
        };
    Status=Prepared;
    INFO("AnnihilationOperator_" << Index <<": " << size << " parts will be computed");
    };
}

BlockNumber FieldOperator::getRightIndex(BlockNumber LeftIndex)
{
	return mapLeftToRightIndex.count(LeftIndex)?mapLeftToRightIndex[LeftIndex]:ERROR_BLOCK_NUMBER;
};

BlockNumber FieldOperator::getLeftIndex(BlockNumber RightIndex)
{
	return mapRightToLeftIndex.count(RightIndex)?mapRightToLeftIndex[RightIndex]:ERROR_BLOCK_NUMBER;
};

QuantumNumbers CreationOperator::mapsTo(QuantumNumbers in) // Require explicit knowledge of QuantumNumbers structure - Not very good
{
  int lz, spin;
  System.getSiteInfo(Index,lz,spin);
  QuantumNumbers q_out;
  if (spin == 1) 
    q_out = QuantumNumbers(in[0] + lz,in[1]+1,in[2]);
  else 
    q_out = QuantumNumbers(in[0] + lz,in[1],in[2]+1);
  return System.checkQuantumNumbers(q_out)?q_out:ERROR_QUANTUM_NUMBERS;
}

BlockNumber CreationOperator::mapsTo(BlockNumber RightIndex)
{
  QuantumNumbers q_right = System.getBlockInfo(RightIndex);	
  QuantumNumbers q_left = (*this).mapsTo(q_right);
  BlockNumber out = System.getBlockNumber(q_left);
  return out;
}

QuantumNumbers AnnihilationOperator::mapsTo(QuantumNumbers in) // Require explicit knowledge of QuantumNumbers structure - Not very good
{
  int lz, spin;
  System.getSiteInfo(Index,lz,spin);
  QuantumNumbers q_out;
  if (spin == 1) 
    q_out = QuantumNumbers(in[0] - lz,in[1]-1,in[2]);
  else 
    q_out = QuantumNumbers(in[0] - lz,in[1],in[2]-1);
  return System.checkQuantumNumbers(q_out)?q_out:ERROR_QUANTUM_NUMBERS;
}

BlockNumber AnnihilationOperator::mapsTo(BlockNumber in)
{
  QuantumNumbers q_in = System.getBlockInfo(in);	
  QuantumNumbers q_out = (*this).mapsTo(q_in);
  BlockNumber out = System.getBlockNumber(q_out);
  return out;
}
