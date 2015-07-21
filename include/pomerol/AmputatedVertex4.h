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


/** \file src/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_VERTEX4_H
#define __INCLUDE_VERTEX4_H

#include"Misc.h"
#include"FourIndexObject.h"
#include"GFContainer.h"
#include"TwoParticleGFContainer.h"

namespace Pomerol{

/** Objects of this class just transforms a two-particle Green's function into
 * an irreducible vertex part or into an amputated irreducible vertex.
 */
class Vertex4 : public FourIndexContainerObject, public Thermal {

    /** A reference to a two-particle Green's function. */
    TwoParticleGFContainer &Chi;
    /** A reference to a Green's function container */
    GFContainer &g;
    /** A reference to a bit classification object */
    const IndexClassification &IndexInfo;
    
    /** Precomputed inverted Green's function matrices calculated at different Matsubara frequencies. */
    std::vector<MatrixType> InvertedGFs;

    /** Amount of computed matsubaras in TwoParticleGF */
    long NumberOfMatsubaras;

    /** A storage for unamputated values */
    std::map<IndexCombination,MatsubaraContainer*> mapUnAmputatedValues;

    /** A vector of all nontrivial combinations to compute */
    std::vector<IndexCombination*> NonTrivialAmputatedCombinations;

    /** A storage for amputated values */
    std::map<IndexCombination,MatsubaraContainer*> mapAmputatedValues;

public:
    /** Constructor.
     * \param[in] IndexInfo A reference to a bit classification object.
     * \param[in] Chi A reference to a two-particle Green's function.
     * \param[in] g1 A reference to a Green's function container.
     */
    Vertex4(const IndexClassification &IndexInfo, TwoParticleGFContainer &Chi, GFContainer &g);

    //============================= UnAmputated methods ==============================//

    /** Do some preparation procedures : prepare storage */
    void prepareUnAmputated();

    /** Compute unamputated values
     */
    void computeUnAmputated();

     /** Returns the value of the unamputated irreducible vertex calculated at given frequencies.
     * \param[in] Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType getUnAmputatedValue(const IndexCombination& in,
                             long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);


    //============================= Amputated methods ==============================//
     /** Do some preparation procedures : calculate inverted GF's, prepare storage */
    void prepareAmputated(std::vector<IndexCombination*>&);

     /** Returns the value of the amputated irreducible vertex calculated at given frequencies.
     * \param[in] Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType getAmputatedValue(const IndexCombination& in,
                             long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    void computeAmputated();

    //==============================================================================//
    /** Returns the value of the irreducible vertex calculated at given frequencies.
     * \param[in] in Requested indices.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(const IndexCombination& in, 
                           long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);


private:
    bool vanishes(const IndexCombination& in); 

    /** Compute unamputated values
     * \param[in] Requested indices.
     */
    void computeUnAmputated(const IndexCombination& in);

    /** Compute amputated values
     * \param[in] Requested indices.
     */
    void computeAmputated(const IndexCombination& in);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_VERTEX4_H
