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


/** \file src/FourIndexObject.h
** \brief A prototype class for all objects, depending on four indices
**
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_FOURINDEXOBJECT_H 
#define __INCLUDE_FOURINDEXOBJECT_H 

#include "Misc.h"
#include "TwoParticleGFPart.h"

namespace Pomerol{

/** This class is a prototype for every object which is defined by a combination of 4 arbitrary Particle Indices
 *  It defines main subobjects, which are used by the 4 index-dependent quantites: IndexCombination and MatsubaraContainer
 */
class FourIndexObject {
public:
    /** A combination of four indices. The notation is ccc^*c^* */
    struct IndexCombination;
    /** A storage of an object over Matsubara frequencies */
    class MatsubaraContainer;
};

/** A prototype container class of Four Indices Objects - stores all values for all possible Particle Indices 
 * The notation is ccc^*c^*
 */ 
class FourIndexContainerObject : public FourIndexObject {
protected:
    /** A set of trivial permutations, which just lead to a change of sign */
    static const Permutation4 TrivialOperatorPermutations[];
public:
    /** Returns the value for a given set of indices and Matsubara numbers (not frequencies themselves)
     * \param[in] In An IndexCombination at which the value should be get.
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
//    ComplexType operator()(const IndexCombination& In, long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
};

/** A structure to handle a combination of 4 Particle Index. The notation is ccc^*c^* */
struct FourIndexObject::IndexCombination
{
    /** Actual indices */
    ParticleIndex Indices[4];
    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& output, const FourIndexObject::IndexCombination& out);
    /** Operator < - comparison method for IndexCombination */
    bool operator< (const FourIndexObject::IndexCombination& rhs) const ;
    /** Operator == */
    bool operator==(const FourIndexObject::IndexCombination& rhs) const ;
    /** Operator != */
    bool operator!=(const FourIndexObject::IndexCombination& rhs) const ;
    /** Constructor
     * \param[in] cindex1 - Index of a 1st operator ( annihilation )
     * \param[in] cindex2 - Index of a 2nd operator ( annihilation )
     * \param[in] cdagindex3 - Index of a 3rd operator ( creation )
     * \param[in] cdagindex4 - Index of a 4th operator ( creation )
     */
    IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_FOURINDEXOBJECT_H 
