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


/** \file Symmetrizer.cpp
**  \brief Implementation of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#include "Symmetrizer.h"

namespace Pomerol { 

//
//Symmetrizer::IndexPermutation
//

Symmetrizer::IndexPermutation::IndexPermutation(const std::vector<ParticleIndex> &in):DynamicIndexCombination(in)
{
    for (ParticleIndex i=0; i<N; ++i) {
        if (in[i]>=N) { ERROR("Indices in IndexPermutation should belong to the interval 0..N-1"); throw (exWrongIndices()); };
        for (ParticleIndex j=i+1; j<N; ++j)
            if (in[i]==in[j]) throw ( exEqualIndices());
        };
    CycleLength=0;
}

const std::vector<ParticleIndex>& Symmetrizer::IndexPermutation::getIndices( unsigned int cycle_number ) const
{
    return Indices;
}

const unsigned int Symmetrizer::IndexPermutation::getCycleLength() const
{
    return CycleLength;
}

const char* Symmetrizer::IndexPermutation::exEqualIndices::what() const throw(){
    return "Cannot have equal indices in the Symmetrizer index combination";
};

//
// Symmetrizer
//

Symmetrizer::Symmetrizer(IndexClassification &IndexInfo):IndexInfo(IndexInfo)
{
}

} // end of namespace Pomerol 

