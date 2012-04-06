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


/** \file Symmetrizer.h
**  \brief Declaration of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_SYMMETRIZER_H
#define __INCLUDE_SYMMETRIZER_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include <set>

namespace Pomerol{

/** This class stores the information about operations, which commute with the Hamiltonian.
 * It tries to find Lattice symmetries, checks for some common symmetries 
 * and also checks given symmetries. */
class Symmetrizer
{
public:
    struct IndexPermutation;
private:
    IndexClassification &IndexInfo;
public:
    Symmetrizer(IndexClassification &IndexInfo);
};

/** A combination of indices to which a permutation commutes with a Hamiltonian. */
struct Symmetrizer::IndexPermutation : protected DynamicIndexCombination
{
private:
    /** This defines how many permutations should be done to form a closed loop. */
    unsigned int CycleLength;
public:
    /** Constructor. 
     * \param[in] @in Equivalent index combination. 
     * \param[in] @length The length of a cycle of permutations.
     */
    IndexPermutation(const std::vector<ParticleIndex> &in);

    /** Returns the permutation for a given number in cycle
     * \param[in] cycle Defines which permutation of indices to return. By default returns the first one.
     */
    const std::vector<ParticleIndex>& getIndices( unsigned int cycle = 1) const;

    /** Returns the length of the permutation cycle. */
    const unsigned int getCycleLength() const;
    /** Exception - equal indices. */
    class exEqualIndices : public std::exception { virtual const char* what() const throw(); };

};

}; // end of namespace Pomerol

#endif //  endif :: #ifndef __INCLUDE_SYMMETRIZER_H
