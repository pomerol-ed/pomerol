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
#include "IndexHamiltonian.h"
#include <set>

namespace Pomerol{

/** This class stores the information about operations, which commute with the Hamiltonian.
 * It tries to find Lattice symmetries, checks for some common symmetries 
 * and also checks given symmetries. */
class Symmetrizer
{
public:
    struct IndexPermutation;
    /** This generates a trivial combination of 0123...N-1 indices. */
    static const DynamicIndexCombination& generateTrivialCombination(ParticleIndex N);
private:
    /** A link to an IndexClassification object. */ 
    const IndexClassification &IndexInfo;
    /** A link to an IndexHamiltonian object. */
    const IndexHamiltonian &Storage;

    bool NSymmetry;
    bool SzSzSymmetry;
    /** A list of equivalent lattice sites permutations. */
    std::list<IndexPermutation*> Permutations;
public:
    Symmetrizer(IndexClassification &IndexInfo, IndexHamiltonian &Storage);
    /** This method finds all lattice permutation operators, that commute with the hamiltonian. */
    //void findLatticeSymmetry();
    /** This method checks the conservation of number of particles. */
    void checkNSymmetry();
    /** This method checks that spin-projection on the z axis is conserved. */
    void checkSzSymmetry();
    /** Returns a list of equivalent permutations. */
    const std::list<IndexPermutation*>& getPermutations() const;
};

/** A combination of indices to which a permutation commutes with a Hamiltonian. 
 * Only the combinations with all different indices can be accepted. Also, since it
 * represents a permutation of indices a check for the irreducibility of the permutation is done, 
 * i.e. it is checked that a permutation can not be split in at least two others. 
 * It is assumed that a trivial identity permutation makes an error. */
struct Symmetrizer::IndexPermutation 
{
private:
    std::vector<DynamicIndexCombination*> Combinations;
    /** This defines how many permutations should be done to form a closed loop. */
    unsigned int CycleLength;
    /** Calculates CycleLength. */
    void calculateCycleLength();
    /** Checks that all elements of permutation are different and belong to 0..N-1 interval. */
    bool checkConsistency(const DynamicIndexCombination &in);
    /** Checks that the permutation can not be splitted. */
    bool checkIrreducibility(const DynamicIndexCombination &in); 
    /** Total amount of indices. */
    const ParticleIndex N;
public:
    /** Constructor. 
     * \param[in] @in Equivalent index combination. 
     */
    IndexPermutation(const DynamicIndexCombination &in);

    /** Returns the permutation for a given number in cycle
     * \param[in] cycle Defines which permutation of indices to return. By default returns the first one.
     */
    const DynamicIndexCombination& getIndices( unsigned int cycle = 1) const;

    /** Returns the length of the permutation cycle. */
    const unsigned int getCycleLength() const;
    /** Exception - equal indices. */
    class exEqualIndices : public std::exception { virtual const char* what() const throw(); };

};

}; // end of namespace Pomerol

#endif //  endif :: #ifndef __INCLUDE_SYMMETRIZER_H
