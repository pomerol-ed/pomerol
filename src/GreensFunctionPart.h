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

    /** Every term is a fraction \f$ \frac{R}{z - P} \f$. */
    struct GreensTerm{
        /** Residue at the pole (\f$ R \f$). */
        ComplexType Residue;
        /** Position of the pole (\f$ P \f$). */
        ComplexType Pole;

        /** Constructor.
         * \param[in] Residue Value of the residue.
         * \param[in] Pole Position of the pole.
         */
        GreensTerm(ComplexType Residue, ComplexType Pole);
        /** Returns a contribution to the Green's function made by this term.
        * \param[in] Frequency Complex frequency \f$ z \f$ to substitute into this term.
        */
        ComplexType operator()(ComplexType Frequency) const;
    };
    /** A stream insertion operator for type GreensTerm.
     * \param[in] out An output stream to insert to.
     * \param[in] Term A term to be inserted.
     */
    friend std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::GreensTerm& Term);

    /** A list of all terms. */
    std::list<GreensTerm> Terms;

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

    /** Iterates over all matrix elements and fills the list of terms. */
    void compute(void);
    /** Returns a sum of all the terms with a substituted Matsubara frequency.
    * \param[in] MatsubaraNum Number of the Matsubara frequency (\f$ \omega_n = \pi*(2*n+1)/\beta \f$).
    */
    ComplexType operator()(long MatsubaraNum) const;
};

std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::GreensTerm& Term);

#endif // endif :: #ifndef __INCLUDE_GREENSFUNCTIONPART_H
