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


/** \file src/GreensFunctionPart.h
** \brief Part of a Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_GREENSFUNCTIONPART_H
#define __INCLUDE_GREENSFUNCTIONPART_H

#include<iomanip>

#include"Misc.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"HamiltonianPart.h"
#include"FieldOperator.h"
#include"DensityMatrixPart.h"

namespace Pomerol{

/** This class represents a part of a Green's function.
 * Every part describes all transitions allowed by selection rules
 * between a given pair of Hamiltonian blocks.
 */
class GreensFunctionPart : public ComputableObject, public Thermal
{
    /** A reference to a part of a Hamiltonian (inner index iterates through it). */
    HamiltonianPart& HpartInner;
    /** A reference to a part of a Hamiltonian (outer index iterates through it). */
    HamiltonianPart& HpartOuter;
    /** A reference to a part of a density matrix (the part corresponding to HpartInner). */
    DensityMatrixPart& DMpartInner;
    /** A reference to a part of a density matrix (the part corresponding to HpartOuter). */
    DensityMatrixPart& DMpartOuter;

    /** A reference to a part of an annihilation operator. */
    AnnihilationOperatorPart& C;
    /** A reference to a part of a creation operator. */
    CreationOperatorPart& CX;

    struct Term;

    /** A stream insertion operator for type GreensTerm.
     * \param[in] out An output stream to insert to.
     * \param[in] Term A term to be inserted.
     */
    friend std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::Term& T);

    /** A list of all terms. */
    std::list<Term> Terms;

    /** A matrix element with magnitude less than this value is treated as zero. */
    static const RealType MatrixElementTolerance = 1e-8;

public:

    /** Constructor.
     * \param[in] C A reference to a part of an annihilation operator.
     * \param[in] CX A reference to a part of a creation operator.
     * \param[in] HpartInner A reference to a part of the Hamiltonian (inner index).
     * \param[in] HpartOuter A reference to a part of the Hamiltonian (outer index).
     * \param[in] DMpartInner A reference to a part of the density matrix (inner index).
     * \param[in] DMpartOuter A reference to a part of the density matrix (outer index).
     */
    GreensFunctionPart(AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                       HamiltonianPart& HpartInner, HamiltonianPart& HpartOuter,
                       DensityMatrixPart& DMpartInner, DensityMatrixPart& DMpartOuter);

    /** Stub prepare() method.*/
    void prepare();

    /** Iterates over all matrix elements and fills the list of terms. */
    void compute(void);
    /** Returns a sum of all the terms with a substituted Matsubara frequency.
    * \param[in] MatsubaraNum Number of the Matsubara frequency (\f$ \omega_n = \pi*(2*n+1)/\beta \f$).
    */
    ComplexType operator()(long MatsubaraNum) const;

    /** Reduces the number of calculated terms 
    * \param[in] Tolerance The tolerance for the terms cutoff.
    * \param[in] ResonantTerms The list of terms.
    */
    static void reduceTerms(const RealType Tolerance, std::list<Term>& Terms);

    /** A difference in energies with magnitude less than this value is treated as zero. */
    static const RealType ReduceResonanceTolerance = 1e-8;//1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
    static const RealType ReduceTolerance = 1e-8;
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

    /** This operator add a term to this one.
    * It does not check the similarity of the terms! 
    * \param[in] AnotherTerm Another term to add to this.
    */
    Term& operator+=(const Term& AnotherTerm);

    /** Returns true if another term is similar to this
     * (sum of the terms is again a correct term).
    */
    bool isSimilarTo(const Term& T) const;
};

std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::Term& T);

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GREENSFUNCTIONPART_H
