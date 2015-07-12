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


/** \file src/GreensFunctionPart.h
** \brief Part of a Green's function.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_GREENSFUNCTIONPART_H
#define __INCLUDE_GREENSFUNCTIONPART_H

#include<iomanip>
#include<cmath>

#include"Misc.h"
#include"StatesClassification.h"
#include"HamiltonianPart.h"
#include"FieldOperator.h"
#include"DensityMatrixPart.h"

namespace Pomerol{

/** This class represents a part of a Green's function.
 * Every part describes all transitions allowed by selection rules
 * between a given pair of Hamiltonian blocks.
 */
class GreensFunctionPart : public Thermal
{
    /** A reference to a part of a Hamiltonian (inner index iterates through it). */
    const HamiltonianPart& HpartInner;
    /** A reference to a part of a Hamiltonian (outer index iterates through it). */
    const HamiltonianPart& HpartOuter;
    /** A reference to a part of a density matrix (the part corresponding to HpartInner). */
    const DensityMatrixPart& DMpartInner;
    /** A reference to a part of a density matrix (the part corresponding to HpartOuter). */
    const DensityMatrixPart& DMpartOuter;

    /** A reference to a part of an annihilation operator. */
    const AnnihilationOperatorPart& C;
    /** A reference to a part of a creation operator. */
    const CreationOperatorPart& CX;

    struct Term;

    /** A stream insertion operator for type GreensTerm.
     * \param[in] out An output stream to insert to.
     * \param[in] Term A term to be inserted.
     */
    friend std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::Term& T);

    /** A list of all terms. */
    std::list<Term> Terms;

    /** A matrix element with magnitude less than this value is treated as zero. */
    const RealType MatrixElementTolerance; // 1e-8;

public:

    /** Constructor.
     * \param[in] C A reference to a part of an annihilation operator.
     * \param[in] CX A reference to a part of a creation operator.
     * \param[in] HpartInner A reference to a part of the Hamiltonian (inner index).
     * \param[in] HpartOuter A reference to a part of the Hamiltonian (outer index).
     * \param[in] DMpartInner A reference to a part of the density matrix (inner index).
     * \param[in] DMpartOuter A reference to a part of the density matrix (outer index).
     */
    GreensFunctionPart(const AnnihilationOperatorPart& C, const CreationOperatorPart& CX, 
                       const HamiltonianPart& HpartInner, const HamiltonianPart& HpartOuter,
                       const DensityMatrixPart& DMpartInner, const DensityMatrixPart& DMpartOuter);

    /** Iterates over all matrix elements and fills the list of terms. */
    void compute(void);

    /** Returns a sum of all the terms with a substituted frequency.
    * \param[in] z Input frequency
    */
    ComplexType operator()(ComplexType z) const;
    /** Returns a sum of all the terms with a substituted Matsubara frequency.
    * \param[in] MatsubaraNum Number of the Matsubara frequency (\f$ \omega_n = \pi*(2*n+1)/\beta \f$).
    */
    ComplexType operator()(long MatsubaraNumber) const;

    /** Returns a sum of all the terms with a substituted imaginary time point.
     * \param[in] tau Imaginary time point.
     */
    ComplexType of_tau(RealType tau) const;

    /** Reduces the number of calculated terms 
    * \param[in] Tolerance The tolerance for the terms cutoff.
    * \param[in] ResonantTerms The list of terms.
    */
    void reduceTerms(const RealType Tolerance, std::list<Term>& Terms);

    /** A difference in energies with magnitude less than this value is treated as zero. */
    const RealType ReduceResonanceTolerance;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
    const RealType ReduceTolerance;
};

/** Every term is a fraction \f$ \frac{R}{z - P} \f$. */
struct GreensFunctionPart::Term {
    /** Residue at the pole (\f$ R \f$). */
    ComplexType Residue;
    /** Position of the pole (\f$ P \f$). */
    RealType Pole;

    /** Constructor.
     * \param[in] Residue Value of the residue.
     * \param[in] Pole Position of the pole.
     */
    Term(ComplexType Residue, RealType Pole);
    /** Returns a contribution to the Green's function made by this term.
     * \param[in] Frequency Complex frequency \f$ z \f$ to substitute into this term.
     */
    ComplexType operator()(ComplexType Frequency) const;

    /** Returns a contribution to the imaginary-time Green's function made by this term.
     * \param[in] tau Imaginary time point.
     * \param[in] beta Inverse temperature.
     */
    ComplexType of_tau(RealType tau, RealType beta) const;

    /** This operator add a term to this one.
    * It does not check the similarity of the terms! 
    * \param[in] AnotherTerm Another term to add to this.
    */
    Term& operator+=(const Term& AnotherTerm);

    /** Returns true if another term is similar to this
     * (sum of the terms is again a correct term).
    */
    bool isSimilarTo(const Term& T, RealType ReduceResonanceTolerance) const;
};

std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::Term& T);

// Inline call operators
inline ComplexType GreensFunctionPart::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing*RealType(2*MatsubaraNumber+1)); }

inline ComplexType GreensFunctionPart::operator()(ComplexType z) const {
    ComplexType G = 0; 
    for(std::list<Term>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm) G += (*pTerm)(z);
    return G;
}

inline ComplexType GreensFunctionPart::of_tau(RealType tau) const {
    ComplexType G = 0; 
    for(std::list<Term>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm) G += pTerm->of_tau(tau,beta);
    return G;
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GREENSFUNCTIONPART_H
