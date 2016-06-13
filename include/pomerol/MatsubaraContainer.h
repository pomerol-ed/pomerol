//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_MATSUBARACONTAINER_H 
#define __INCLUDE_MATSUBARACONTAINER_H 

#include "Misc.h"

namespace Pomerol{

class MatsubaraContainer 
{
    /** This is a i\frac{\pi}{\beta} - an interval between 2 adjacent matsubaras. 
     * It is more convenient to store Inverse Temperature with a spacing, rather than to invert it afterwards 
     */
    ComplexType MatsubaraSpacing;

    /* An amount of positive Matsubaras on which to store the values. The range of values [-MatsubaraSpacing;MatsubaraSpacing-1] will be available */
    int FermionicMin_;
    int FermionicMax_;
    int BosonicMin_;
    int BosonicMax_;

    /* The storage - an array of Matrices for nu,nu' space, stored as a vector which is dependent on a bosonic frequency index */
    std::vector<ComplexMatrixType> Data;
    /* A first non-vanishing fermionic index - used for correspondence between number of element in a matrix and real Matsubara numbers */
    std::vector<long> FermionicFirstIndex;
public:
    /** Constructor
     * \param[in] beta An inverse temperature.
     */
    MatsubaraContainer(RealType beta);

    //void prepare(int NBosonic, int NFermionic);
    /** Allocate memory for a storage
     * \param[in] NumberOfMatsubaras An amount of positive Matsubara frequencies that will be held in a MatsubaraContainer.
     */
    void prepare(int BosonicMin, int BosonicMax, int FermionicMin, int FermionicMax);

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
    //void set(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3, ComplexType &Value);


    /** Fill container from a list of terms 
     * \param[in] NonResonantTerms A list of NonResonant Terms.
     * \param[in] ResonantTerms A list of Resonant Terms.
     * \param[in] Permutation A permutation of input Matsubara frequencies
     */
    //void fill(const std::vector<TwoParticleGFPart::NonResonantTerm>& NonResonantTerms, const std::vector<TwoParticleGFPart::ResonantTerm>& ResonantTerms, Permutation3 Permutation);

    /** Operator+= : adds to a current MatsubaraContainer another one
     * \param[in] rhs Right hand side of the equation Matsubara Container to add.
     */
    MatsubaraContainer& operator+=(const MatsubaraContainer& rhs);

    /** Empty memory */
    void clear();

    int NBosonic() const { return BosonicMax_ - BosonicMin_ + 1; }
    int NFermionic() const { return FermionicMax_ - FermionicMin_;}

    friend class TwoParticleGF;
    friend class TwoParticleGFPart;
};

inline
ComplexType MatsubaraContainer::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
// w_1 = w                w = w_1
// w_2 = w' + W           W = w_3 - w_1
// w_3 = w + W            w'= w_2 - W
{
    int nu1 = MatsubaraNumber1;
    int W = MatsubaraNumber3 - MatsubaraNumber1;
    int nu2 = MatsubaraNumber2 - W;

    int WIndex = W - BosonicMin_;
    int nu1Index = nu1 - FermionicMin_;
    int nu2Index = nu2 - FermionicMin_;

    //cout << "Bosonic index : " << RealBosonicIndex - 2*NumberOfMatsubaras<< " shift : " << FermionicFirstIndex[RealBosonicIndex] << endl;
    if (WIndex >= 0 && WIndex <= BosonicMax_ - BosonicMin_ 
        && nu1Index >= 0 && nu1 <= FermionicMax_   
        && nu2Index >= 0 && nu2 <= FermionicMax_)
        return Data[WIndex](nu1Index,nu2Index);
    else {
        ERROR("Warning! Matsubara numbers (" << MatsubaraNumber1 << "," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+ MatsubaraNumber2 - MatsubaraNumber3 << ") of FourIndexObject is out of range, returning 0");
        return ComplexType (0.0,0.0);
    };
};

/*
inline
void MatsubaraContainer::set(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3, ComplexType &Value)
{
    unsigned int RealBosonicIndex = MatsubaraNumber1 + MatsubaraNumber2 + 2*NumberOfMatsubaras;
    int nuIndex = MatsubaraNumber1-FermionicFirstIndex[RealBosonicIndex];
    int nu1Index= RealBosonicIndex-2*NumberOfMatsubaras-MatsubaraNumber3-FermionicFirstIndex[RealBosonicIndex];
    if (nuIndex >= 0 && nuIndex < Data[RealBosonicIndex].rows() && nu1Index >= 0 && nu1Index < Data[RealBosonicIndex].rows() )
        Data[RealBosonicIndex](nuIndex,nu1Index)=Value;
    else ERROR("Warning! Tried assigning to wrong Matsubara numbers (" << MatsubaraNumber1 << "," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+ MatsubaraNumber2 - MatsubaraNumber3 <<"). Value left unassigned");
}
*/



} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_MATSUBARACONTAINER_H
