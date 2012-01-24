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


/**
* A subclass of FourIndexObject - it stores values of a 4 FourIndexObject  over Matsubara frequencies. 
* The data is stored in a 1 bosonic and 2 fermionic degrees flavour : 
* (OMEGA,nu,nu'), 
* where OMEGA=w1+w2 - bosonic frequency, nu=w1, nu'=w4
*/
class FourIndexObject::MatsubaraContainer 
{
    /** This is a i\frac{\pi}{\beta} - an interval between 2 adjacent matsubaras. 
     * It is more convenient to store Inverse Temperature with a spacing, rather than to invert it afterwards 
     */
    ComplexType MatsubaraSpacing;

    /* An amount of positive Matsubaras on which to store the values. The range of values [-MatsubaraSpacing;MatsubaraSpacing-1] will be available */
    long NumberOfMatsubaras;

    /* The storage - an array of Matrices for nu,nu' space, stored as a vector which is dependent on a bosonic frequency index */
    std::vector<MatrixType> Data;
    /* A first non-vanishing fermionic index - used for correspondence between number of element in a matrix and real Matsubara numbers */
    std::vector<long> FermionicFirstIndex;
public:
    /** Constructor
     * \param[in] beta An inverse temperature.
     */
    MatsubaraContainer(RealType beta);

    /** Allocate memory for a storage
     * \param[in] NumberOfMatsubaras An amount of positive Matsubara frequencies that will be held in a MatsubaraContainer.
     */
    void prepare(long NumberOfMatsubaras);

    /** Returns the value for a given Matsubara numbers (not frequencies themselves)
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /** Sets the value for a given Matsubara numbers (not frequencies themselves)
     * \param[in] MatsubaraNumber1 An index of the 1st Matsubara frequency.
     * \param[in] MatsubaraNumber2 An index of the 2nd Matsubara frequency.
     * \param[in] MatsubaraNumber3 An index of the 3rd Matsubara frequency.
     */
    void set(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3, ComplexType &Value);


    /** Fill container from a list of terms 
     * \param[in] NonResonantTerms A list of NonResonant Terms.
     * \param[in] ResonantTerms A list of Resonant Terms.
     * \param[in] Permutation A permutation of input Matsubara frequencies
     */
    void fill(const std::list<TwoParticleGFPart::NonResonantTerm>& NonResonantTerms, const std::list<TwoParticleGFPart::ResonantTerm>& ResonantTerms, Permutation3 Permutation);

    /** Operator+= : adds to a current MatsubaraContainer another one
     * \param[in] rhs Right hand side of the equation Matsubara Container to add.
     */
    MatsubaraContainer& operator+=(const MatsubaraContainer& rhs);

    /** Empty memory */
    void clear();
};

inline
ComplexType FourIndexObject::MatsubaraContainer::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
// {OMEGA,nu,nu' : OMEGA=w1+w2, nu=w1, nu'=w4=OMEGA-w3
{
    unsigned int RealBosonicIndex = MatsubaraNumber1 + MatsubaraNumber2 + 2*NumberOfMatsubaras;
    int nuIndex = MatsubaraNumber1-FermionicFirstIndex[RealBosonicIndex];
    int nu1Index= RealBosonicIndex-2*NumberOfMatsubaras-MatsubaraNumber3-FermionicFirstIndex[RealBosonicIndex];
    //cout << "Bosonic index : " << RealBosonicIndex - 2*NumberOfMatsubaras<< " shift : " << FermionicFirstIndex[RealBosonicIndex] << endl;
    if (nuIndex >= 0 && nuIndex < Data[RealBosonicIndex].rows() && nu1Index >= 0 && nu1Index < Data[RealBosonicIndex].rows() )
        return Data[RealBosonicIndex](nuIndex,nu1Index);
    else {
        ERROR("Warning! Matsubara numbers (" << MatsubaraNumber1 << "," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+ MatsubaraNumber2 - MatsubaraNumber3 << ") of FourIndexObject is out of range, returning 0");
        return ComplexType (0.0,0.0);
    };
};

inline
void FourIndexObject::MatsubaraContainer::set(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3, ComplexType &Value)
{
    unsigned int RealBosonicIndex = MatsubaraNumber1 + MatsubaraNumber2 + 2*NumberOfMatsubaras;
    int nuIndex = MatsubaraNumber1-FermionicFirstIndex[RealBosonicIndex];
    int nu1Index= RealBosonicIndex-2*NumberOfMatsubaras-MatsubaraNumber3-FermionicFirstIndex[RealBosonicIndex];
    if (nuIndex >= 0 && nuIndex < Data[RealBosonicIndex].rows() && nu1Index >= 0 && nu1Index < Data[RealBosonicIndex].rows() )
        Data[RealBosonicIndex](nuIndex,nu1Index)=Value;
    else ERROR("Warning! Tried assigning to wrong Matsubara numbers (" << MatsubaraNumber1 << "," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+ MatsubaraNumber2 - MatsubaraNumber3 <<"). Value left unassigned");
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_FOURINDEXOBJECT_H 
