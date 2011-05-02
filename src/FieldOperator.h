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


#ifndef __INCLUDE_FIELDOPERATOR_H
#define __INCLUDE_FIELDOPERATOR_H

#include"Misc.h"
#include"ComputableObject.h"
#include"output.h"
#include"StatesClassification.h"
#include"Hamiltonian.h" 
#include"FieldOperatorPart.h"

typedef std::pair<BlockNumber,BlockNumber> BlockMapping;

class FieldOperator : public ComputableObject
{
protected:
    IndexClassification &IndexInfo;
    StatesClassification &System;
    Hamiltonian &H;

    ParticleIndex Index;
    std::vector<FieldOperatorPart*> Data;
    std::map<unsigned int,BlockNumber> mapPartsFromRight;        // A map from non-zero parts to their BlockNumber indices
    std::map<unsigned int,BlockNumber> mapPartsFromLeft;        // A map from output index to input index, hence there is a unique transform
    std::map<BlockNumber,BlockNumber> mapRightToLeftIndex;        // A map from output index to input index, hence there is a unique transform
    std::map<BlockNumber,BlockNumber> mapLeftToRightIndex;        // A map from output index to input index, hence there is a unique transform
    std::list<BlockMapping> LeftRightIndices;

    virtual BlockNumber mapsTo(BlockNumber RightIndex)=0;
    virtual QuantumNumbers mapsTo(QuantumNumbers in)=0;

public:
    FieldOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index);

    FieldOperatorPart& getPartFromLeftIndex(BlockNumber in);
    FieldOperatorPart& getPartFromLeftIndex(QuantumNumbers in);
    FieldOperatorPart& getPartFromRightIndex(BlockNumber out);
    FieldOperatorPart& getPartFromRightIndex(QuantumNumbers out);
    BlockNumber getLeftIndex(BlockNumber RightIndex);
    BlockNumber getRightIndex(BlockNumber LeftIndex);
    std::list<BlockMapping>& getNonTrivialIndices();

    virtual void prepare()=0;
    
    void compute();
    void dump();
    void print_to_screen();
    unsigned short getIndex() const;
};

class CreationOperator;
class AnnihilationOperator;

class CreationOperator : public FieldOperator
{
    friend class AnnihilationOperator;
    BlockNumber mapsTo(BlockNumber RightIndex);
    QuantumNumbers mapsTo(QuantumNumbers in);
public:
    AnnihilationOperator& transpose();
        void prepare();
    
    CreationOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index);
};

class AnnihilationOperator : public FieldOperator
{
    friend class CreationOperator;
    BlockNumber mapsTo(BlockNumber RightIndex);
    QuantumNumbers mapsTo(QuantumNumbers in);
public:
    CreationOperator& transpose();
        void prepare();
    
    AnnihilationOperator(IndexClassification &IndexInfo, StatesClassification &System, Hamiltonian &H, ParticleIndex Index);
};

#endif // endif :: #ifdef __INCLUDE_FIELDOPERATOR_H
