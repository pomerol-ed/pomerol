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


/** \file src/FieldOperatorContainer.cpp
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "FieldOperatorContainer.h"

FieldOperatorContainer::FieldOperatorContainer(StatesClassification &S, Hamiltonian &H, IndexClassification &IndexInfo):S(S),H(H),IndexInfo(IndexInfo)
{
};

CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapCreationOperators.count(in)==0){
        CreationOperator *CX = new CreationOperator(S,H,in);
        INFO("FieldOperatorContainer: Making Creation Operator_"<<in);
        CX->prepare();
        CX->compute();
        mapCreationOperators[in] = CX;
        };
    //else INFO("FieldOperatorContainer: Using already computed Creation Operator_"<< in);
    return *mapCreationOperators[in];
    }
else assert(0);
}

AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in)
{
if (IndexInfo.checkIndex(in)){
    if (mapAnnihilationOperators.count(in)==0){
        AnnihilationOperator *C = new AnnihilationOperator(S,H,in);
        INFO("FieldOperatorContainer: Making Annihilation Operator_"<<in);
        C->prepare();
        C->compute();
        mapAnnihilationOperators[in] = C;
        };
   // else INFO("FieldOperatorContainer: Using already computed Annihilation Operator_"<< in);
    return *mapAnnihilationOperators[in];
    }
else assert(0);
}

