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

#include "pomerol/FieldOperatorContainer.h"

namespace Pomerol{

FieldOperatorContainer::FieldOperatorContainer(IndexClassification &IndexInfo, StatesClassification &S, const Hamiltonian &H, bool use_transpose) : 
    IndexInfo(IndexInfo), S(S), H(H), use_transpose(use_transpose)
{}

void FieldOperatorContainer::prepareAll(std::set<ParticleIndex> in)
{
    if (in.size() == 0) for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i) in.insert(i);
    for (std::set<ParticleIndex>::const_iterator it = in.begin(); it!=in.end(); it++)
        {
            ParticleIndex i = *it;
            CreationOperator *CX = new CreationOperator(IndexInfo, S,H,i);
            CX->prepare();
            mapCreationOperators[i] = CX;
            AnnihilationOperator *C = new AnnihilationOperator(IndexInfo, S,H,i);
            C->prepare();
            mapAnnihilationOperators[i] = C;
        }
}

void FieldOperatorContainer::computeAll()
{
    for (std::map <ParticleIndex, CreationOperator*>::iterator cdag_it = mapCreationOperators.begin(); cdag_it != mapCreationOperators.end(); ++cdag_it) {
        CreationOperator &cdag = *(cdag_it->second);
        cdag.compute();
        AnnihilationOperator &c = *mapAnnihilationOperators[cdag_it->first];

        FieldOperator::BlocksBimap cdag_block_map = cdag.getBlockMapping();
        // hack - copy transpose matrices into c
        for (FieldOperator::BlocksBimap::right_const_iterator cdag_map_it=cdag_block_map.right.begin(); cdag_map_it!=cdag_block_map.right.end(); cdag_map_it++) {
                c.getPartFromRightIndex(cdag_map_it->second).elementsRowMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getColMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).elementsColMajor = cdag.getPartFromRightIndex(cdag_map_it->first).getRowMajorValue().adjoint();
                c.getPartFromRightIndex(cdag_map_it->second).Status = ComputableObject::Computed;
                c.Status = ComputableObject::Computed;
            };
        };

// original
    //for (auto c : mapAnnihilationOperators) c.second->compute();
}

const CreationOperator& FieldOperatorContainer::getCreationOperator(ParticleIndex in) const
{
    if (IndexInfo.checkIndex(in)){
        return *mapCreationOperators[in];
        }
    else
        throw (std::logic_error("No creation operator found."));
}

const AnnihilationOperator& FieldOperatorContainer::getAnnihilationOperator(ParticleIndex in) const
{
    if (IndexInfo.checkIndex(in)){
        return *mapAnnihilationOperators[in];
        }
    else
        throw (std::logic_error("No annihilation operator found."));
}

} // end of namespace Pomerol
