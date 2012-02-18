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


#ifndef __INCLUDE_FIELDOPERATORPART_H
#define __INCLUDE_FIELDOPERATORPART_H
#include"Misc.h"
#include"StatesClassification.h"
#include"Hamiltonian.h"

namespace Pomerol{

class FieldOperatorPart {
protected:

    const IndexClassification &IndexInfo;
    const StatesClassification &S;
    const HamiltonianPart &HFrom;
    const HamiltonianPart &HTo;

    ParticleIndex PIndex;
    RowMajorMatrixType elementsRowMajor;
    ColMajorMatrixType elementsColMajor;
    // basic functions
#warning The following 3 functions need more descriptive names... and bodies.
    virtual QuantumState retK(QuantumState L) const = 0;  
    virtual int mFunc(QuantumState state1, QuantumState state2, ParticleIndex i) const = 0;   //checks matrix element of an operator between state1 and state2
    virtual bool checkL(QuantumState L) const = 0; //checks state L to be appropriate as a result of a creation/destruction operator

    static const RealType MatrixElementTolerance = 1e-8;

public:

    FieldOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);

    void compute();
    // DEPRECATED
    void print_to_screen();                        //print to screen matrices UXCU UXCXU

    const RowMajorMatrixType& getRowMajorValue(void) const;
    const ColMajorMatrixType& getColMajorValue(void) const;
    BlockNumber getLeftIndex(void) const;
    BlockNumber getRightIndex(void) const;
};

class AnnihilationOperatorPart;
class CreationOperatorPart;

class AnnihilationOperatorPart : public FieldOperatorPart
{ 
    QuantumState retK(QuantumState L) const;
    int mFunc(QuantumState state1, QuantumState state2, ParticleIndex PIndex) const;
    bool checkL(QuantumState L) const;
    friend class CreationOperatorPart;

public :
    AnnihilationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);
    const CreationOperatorPart& transpose(void) const;
};

class CreationOperatorPart : public FieldOperatorPart
{
    QuantumState retK(QuantumState L) const;
    int mFunc(QuantumState state1, QuantumState state2, ParticleIndex PIndex) const;
    bool checkL(QuantumState L) const;
    friend class AnnihilationOperatorPart;

public :
    CreationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);
    const AnnihilationOperatorPart& transpose(void) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H
