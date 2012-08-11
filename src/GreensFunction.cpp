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


/** \file src/GreensFunction.cpp
** \brief Thermal Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "GreensFunction.h"

namespace Pomerol{

GreensFunction::GreensFunction(const StatesClassification& S, const Hamiltonian& H, 
                               const AnnihilationOperator& C, const CreationOperator& CX,
                               const DensityMatrix& DM) :
    Thermal(DM), ComputableObject(Constructed), S(S), H(H), C(C), CX(CX), DM(DM), Vanishing(true),
    parts(0), Storage(this)
{
}

GreensFunction::~GreensFunction()
{
    for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

void GreensFunction::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of C and CX.
    const std::list<BlockMapping> &CNontrivialBlocks = C.getNonTrivialIndices();
    const std::list<BlockMapping> &CXNontrivialBlocks = CX.getNonTrivialIndices();

    std::list<BlockMapping>::const_iterator Citer = CNontrivialBlocks.begin();
    std::list<BlockMapping>::const_iterator CXiter = CXNontrivialBlocks.begin();

    while(Citer != CNontrivialBlocks.end() && CXiter != CXNontrivialBlocks.end()){
        // <Cleft|C|Cright><CXleft|CX|CXright>
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->first;
        BlockNumber CXright = CXiter->second;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(Cleft == CXright && Cright == CXleft){
        parts.push_back(new GreensFunctionPart(
                              (AnnihilationOperatorPart&)C.getPartFromLeftIndex(Cleft),
                              (CreationOperatorPart&)CX.getPartFromRightIndex(CXright),
                              H.getPart(Cright), H.getPart(Cleft),
                              DM.getPart(Cright), DM.getPart(Cleft)));
        }

        unsigned long CleftInt = Cleft;
        unsigned long CXrightInt = CXright;

        if(CleftInt <= CXrightInt) Citer++;
        if(CleftInt >= CXrightInt) CXiter++;
    }
    if (parts.size() > 0) Vanishing = false;

    Status = Prepared;
}

void GreensFunction::compute(long NumberOfMatsubaras)
{
    if(Status>=Computed) return;
    if(Status<Prepared) prepare();

    if(Status<Computed){
        for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
            (*iter)->compute();
    }
    if(NumberOfMatsubaras != Storage.getNumberOfMatsubaras())
        Storage.fill(NumberOfMatsubaras);
    Status = Computed;
}

inline
ComplexType GreensFunction::value(long int MatsubaraNum) const
{
    ComplexType Value = 0;
    for(std::list<GreensFunctionPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
        Value += (**iter)(MatsubaraNum);
    return Value;
}

ComplexType GreensFunction::operator()(long MatsubaraNum) const
{
    if(Vanishing)
        return 0;
    else
        return Storage(MatsubaraNum);
}

unsigned short GreensFunction::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C.getIndex();
        case 1: return CX.getIndex();
        default: assert(0);
    }
}

bool GreensFunction::isVanishing(void) const
{
    return Vanishing;
}

} // end of namespace Pomerol
