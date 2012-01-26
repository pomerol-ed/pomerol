//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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

namespace Pomerol{
FieldOperator::FieldOperator(IndexClassification &IndexInfo, StatesClassification &System, const Hamiltonian &H, ParticleIndex Index) :
    IndexInfo(IndexInfo), System(System), H(H), Index(Index)
{}

CreationOperator::CreationOperator(IndexClassification &IndexInfo, StatesClassification &System, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,System,H,Index)
{}

AnnihilationOperator::AnnihilationOperator(IndexClassification &IndexInfo, StatesClassification &System, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,System,H,Index)
{}

const std::list<BlockMapping>& FieldOperator::getNonTrivialIndices() const
{
    return LeftRightIndices;
}
 
const FieldOperatorPart& FieldOperator::getPartFromRightIndex(BlockNumber in) const
{
    return *parts[mapPartsFromRight.find(in)->second];
}

const FieldOperatorPart& FieldOperator::getPartFromRightIndex(const QuantumNumbers& in) const
{
    return *parts[mapPartsFromRight.find(System.getBlockNumber(in))->second];
}

const FieldOperatorPart& FieldOperator::getPartFromLeftIndex(BlockNumber in) const
{
    return *parts[mapPartsFromLeft.find(in)->second];
}

const FieldOperatorPart& FieldOperator::getPartFromLeftIndex(const QuantumNumbers& in) const
{
    return *parts[mapPartsFromLeft.find(System.getBlockNumber(in))->second];
}

void FieldOperator::compute(void)
{
    size_t Size = parts.size();
    INFO_NONEWLINE("FieldOperator_" << Index << ", computing: ")
    for (size_t BlockIn = 0; BlockIn < Size; BlockIn++){
        INFO_NONEWLINE( (int) ((1.0*BlockIn/Size) * 100 ) << "  " << std::flush);
        parts[BlockIn]->compute();
    };
    INFO("");
}

ParticleIndex FieldOperator::getIndex(void) const
{
    return Index;
}

void CreationOperator::prepare(void)
{
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0; RightIndex<System.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new CreationOperatorPart(IndexInfo, System,
                                    H.getPart(RightIndex),H.getPart(LeftIndex),Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightIndices.push_back(BlockMapping(LeftIndex,RightIndex));
            mapRightToLeftIndex[RightIndex]=LeftIndex;
            mapLeftToRightIndex[LeftIndex]=RightIndex;
            Size++;
        }
    }
    INFO("CreationOperator_" << Index <<": " << Size << " parts will be computed");
}

void AnnihilationOperator::prepare()
{
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0;RightIndex<System.NumberOfBlocks();RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new AnnihilationOperatorPart(IndexInfo, System,
                                    H.getPart(RightIndex),H.getPart(LeftIndex), Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            mapRightToLeftIndex[RightIndex]=LeftIndex;
            mapLeftToRightIndex[LeftIndex]=RightIndex;
            LeftRightIndices.push_back(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("AnnihilationOperator_" << Index <<": " << Size << " parts will be computed");
}

BlockNumber FieldOperator::getRightIndex(BlockNumber LeftIndex) const
{
    return mapLeftToRightIndex.count(LeftIndex) ?
        mapLeftToRightIndex.find(LeftIndex)->second : ERROR_BLOCK_NUMBER;
}

BlockNumber FieldOperator::getLeftIndex(BlockNumber RightIndex) const
{
    return mapRightToLeftIndex.count(RightIndex) ? 
        mapRightToLeftIndex.find(RightIndex)->second : ERROR_BLOCK_NUMBER;
}

QuantumNumbers CreationOperator::mapsTo(const QuantumNumbers& in) const // Require explicit knowledge of QuantumNumbers structure - Not very good
{
    int lz, spin;
    System.getSiteInfo(Index,lz,spin);
    QuantumNumbers q_out;
    if (spin == 1) 
        q_out = QuantumNumbers(in[0] + lz,in[1]+1,in[2]);
    else 
        q_out = QuantumNumbers(in[0] + lz,in[1],in[2]+1);
    return System.checkQuantumNumbers(q_out) ? q_out : ERROR_QUANTUM_NUMBERS;
}

BlockNumber CreationOperator::mapsTo(BlockNumber RightIndex) const
{
    QuantumNumbers q_right = System.getBlockInfo(RightIndex);	
    QuantumNumbers q_left = mapsTo(q_right);
    BlockNumber out = System.getBlockNumber(q_left);
    return out;
}

QuantumNumbers AnnihilationOperator::mapsTo(const QuantumNumbers& in) const // Require explicit knowledge of QuantumNumbers structure - Not very good
{
    int lz, spin;
    System.getSiteInfo(Index,lz,spin);
    QuantumNumbers q_out;
    if (spin == 1) 
        q_out = QuantumNumbers(in[0] - lz,in[1]-1,in[2]);
    else 
        q_out = QuantumNumbers(in[0] - lz,in[1],in[2]-1);
    return System.checkQuantumNumbers(q_out) ? q_out : ERROR_QUANTUM_NUMBERS;
}

BlockNumber AnnihilationOperator::mapsTo(BlockNumber in) const
{
    QuantumNumbers q_in = System.getBlockInfo(in);	
    QuantumNumbers q_out = mapsTo(q_in);
    BlockNumber out = System.getBlockNumber(q_out);
    return out;
}

} // end of namespace Pomerol
