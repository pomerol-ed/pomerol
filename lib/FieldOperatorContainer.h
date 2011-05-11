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


/** \file src/FieldOperatorContainer.h
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_FIELDOPERATORCONTAINER_H
#define __INCLUDE_FIELDOPERATORCONTAINER_H

#include"Misc.h"
#include"ComputableObject.h"
#include"FieldOperator.h"
#include"Hamiltonian.h"
#include"StatesClassification.h"

/** This class represents a container to store and retrieve FieldOperators ( CreationOperator or AnnihilationOperator 
 * rotated to eigenvector basis of Hamiltonian H ) for a given Index.
 * If no field operator is yet initialized then calculation of the field operator is done.
 */
class FieldOperatorContainer
{
private:
    /** A reference to a states classification object. */
    StatesClassification &S;
    /** A reference to a Hamiltonian. */
    Hamiltonian &H;
    /** A reference to a IndexClassification object in order to check the input indices. */
    IndexClassification &IndexInfo;
    /** A map which gives a link to the CreationOperator for a given index */
    std::map <ParticleIndex, CreationOperator*>     mapCreationOperators;
    /** A map which gives a link to the AnnihilationOperator for a given index */
    std::map <ParticleIndex, AnnihilationOperator*> mapAnnihilationOperators;
public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] IndexInfo A reference to a IndexClassification
     */
    FieldOperatorContainer(StatesClassification &S, 
        Hamiltonian &H, IndexClassification &IndexInfo);

    /** Returns the CreationOperator for a given Index */
    CreationOperator& getCreationOperator(ParticleIndex in);
    /** Returns the AnnihilationOperator for a given Index */
    AnnihilationOperator& getAnnihilationOperator(ParticleIndex in);
};

#endif // endif :: #ifndef __INCLUDE_FIELDOPERATORCONTAINER_H
