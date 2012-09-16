//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "FieldOperatorContainer.h"

namespace Pomerol{

FieldOperatorContainer::FieldOperatorContainer(IndexClassification &IndexInfo, StatesClassification &S, const Hamiltonian &H, bool use_transpose) : 
    ComputableObject(Constructed), IndexInfo(IndexInfo), S(S), H(H), use_transpose(use_transpose)
{}

void FieldOperatorContainer::prepare()
{
    if ( Status >= Prepared ) return;
    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i)
        {
            CreationOperator *CX = new CreationOperator(IndexInfo, S,H,i);
            CX->prepare();
            mapCreationOperators[i] = CX;
            //if (!use_transpose) {
            AnnihilationOperator *C = new AnnihilationOperator(IndexInfo, S,H,i);
            C->prepare();
            mapAnnihilationOperators[i] = C;
           //     }
        }

    Status = Prepared;
}

const CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in) const
{
    if (Status<Prepared) { ERROR("GFContainer needs to be prepared."); throw (exStatusMismatch()); }
    if (IndexInfo.checkIndex(in)){
        mapCreationOperators[in]->compute();
        return *mapCreationOperators[in];
        }
    else
        assert(0);
}

const AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in) const
{
    if (Status<Prepared) { ERROR("GFContainer needs to be prepared."); throw (exStatusMismatch()); }
    if (IndexInfo.checkIndex(in)){
        mapAnnihilationOperators[in]->compute();
        return *mapAnnihilationOperators[in];
        }
    else
        assert(0);
}

} // end of namespace Pomerol
