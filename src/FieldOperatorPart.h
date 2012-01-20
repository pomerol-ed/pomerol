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


#ifndef __INCLUDE_FIELDOPERATORPART_H
#define __INCLUDE_FIELDOPERATORPART_H
#include"Misc.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"Hamiltonian.h"

namespace Pomerol{

class FieldOperatorPart : public ComputableObject {
protected:

    IndexClassification &IndexInfo;
    StatesClassification &S;
    HamiltonianPart &h_from;
    HamiltonianPart &h_to;

    ParticleIndex i;
    RowMajorMatrixType elementsRowMajor;    
    ColMajorMatrixType elementsColMajor;    
    // basic functions
    virtual QuantumState retK(QuantumState L)=0;  
    virtual int mFunc(QuantumState state1, QuantumState state2, ParticleIndex i)=0;   //checks matrix element of an operator between state1 and state2
    virtual bool checkL(QuantumState L)=0;  //checks state L to be appropriate as a result of a creation/destruction operator
    
    static const RealType MatrixElementTolerance = 1e-8;
   
public:
  
    FieldOperatorPart(IndexClassification &IndexInfo, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, ParticleIndex i);

    void prepare();
    void compute();
    void dump();
    void print_to_screen();                        //print to screen matrices UXCU UXCXU

    RowMajorMatrixType& getRowMajorValue();
    ColMajorMatrixType& getColMajorValue();
    BlockNumber getLeftIndex();
    BlockNumber getRightIndex();
};

class AnnihilationOperatorPart;
class CreationOperatorPart;

class AnnihilationOperatorPart : public FieldOperatorPart
{ 
    QuantumState retK(QuantumState L);    
    int mFunc(QuantumState state1, QuantumState state2, ParticleIndex i);
    bool checkL(QuantumState L);
    friend class CreationOperatorPart;

public :
    AnnihilationOperatorPart(IndexClassification &IndexInfo, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, ParticleIndex i);
    CreationOperatorPart& transpose();
};

class CreationOperatorPart : public FieldOperatorPart
{
    QuantumState retK(QuantumState L);    
    int mFunc(QuantumState state1, QuantumState state2, ParticleIndex i);    
    bool checkL(QuantumState L);
    friend class AnnihilationOperatorPart;
  
public :
    CreationOperatorPart(IndexClassification &IndexInfo, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, ParticleIndex i);
    AnnihilationOperatorPart& transpose();
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H
