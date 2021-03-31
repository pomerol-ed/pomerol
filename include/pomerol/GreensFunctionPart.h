/** \file include/pomerol/GreensFunctionPart.h
** \brief Part of a Green's function for a given set of quantum numbers.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_GREENSFUNCTIONPART_H
#define __INCLUDE_GREENSFUNCTIONPART_H

#include<iomanip>
#include<cmath>

#include <mpi.h>

#include"Misc.h"
#include"StatesClassification.h"
#include"HamiltonianPart.h"
#include"FieldOperator.h"
#include"DensityMatrixPart.h"
#include"TermList.h"

namespace Pomerol{

/** This class represents a part of a Green's function.
 * Every part describes all transitions allowed by selection rules
 * between a given pair of Hamiltonian blocks.
 */
template<bool Complex = false>
class GreensFunctionPart : public Thermal
{
    /** A reference to a part of a Hamiltonian (inner index iterates through it). */
    const HamiltonianPart<Complex>& HpartInner;
    /** A reference to a part of a Hamiltonian (outer index iterates through it). */
    const HamiltonianPart<Complex>& HpartOuter;
    /** A reference to a part of a density matrix (the part corresponding to HpartInner). */
    const DensityMatrixPart<Complex>& DMpartInner;
    /** A reference to a part of a density matrix (the part corresponding to HpartOuter). */
    const DensityMatrixPart<Complex>& DMpartOuter;

    /** A reference to a part of an annihilation operator. */
    const AnnihilationOperatorPart<Complex>& C;
    /** A reference to a part of a creation operator. */
    const CreationOperatorPart<Complex>& CX;

    /** Every term is a fraction \f$ \frac{R}{z - P} \f$. */
    struct Term {
        /** Residue at the pole (\f$ R \f$). */
        ComplexType Residue;
        /** Position of the pole (\f$ P \f$). */
        RealType Pole;

        /** Comparator object for terms */
        struct Compare {
            const double Tolerance;
            Compare(double Tolerance) : Tolerance(Tolerance) {}
            bool operator()(Term const& t1, Term const& t2) const {
                return t2.Pole - t1.Pole >= Tolerance;
            }
        };

        /** Does term have a negligible residue? */
        struct IsNegligible {
            double Tolerance;
            IsNegligible(double Tolerance) : Tolerance(Tolerance) {}
            bool operator()(Term const& t, size_t ToleranceDivisor) const {
                return std::abs(t.Residue) < Tolerance / ToleranceDivisor;
            }
            void broadcast(const MPI_Comm &comm, int root) {
                MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm);
            }
        };

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
        ComplexType operator()(RealType tau, RealType beta) const;

        /** This operator add a term to this one.
        * It does not check the similarity of the terms!
        * \param[in] AnotherTerm Another term to add to this.
        */
        Term& operator+=(const Term& AnotherTerm);
    };
    /** A stream insertion operator for type Term.
     * \param[in] out An output stream to insert to.
     * \param[in] Term A term to be inserted.
     */
    friend std::ostream& operator<< (std::ostream& out, const Term& T) {
        out << T.Residue << "/(z - " << T.Pole << ")";
        return out;
    }

    /** A list of all terms. */
    TermList<Term> Terms;

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
    GreensFunctionPart(const AnnihilationOperatorPart<Complex>& C,
                       const CreationOperatorPart<Complex>& CX,
                       const HamiltonianPart<Complex>& HpartInner,
                       const HamiltonianPart<Complex>& HpartOuter,
                       const DensityMatrixPart<Complex>& DMpartInner,
                       const DensityMatrixPart<Complex>& DMpartOuter);

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

    /** A difference in energies with magnitude less than this value is treated as zero. */
    const RealType ReduceResonanceTolerance;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
    const RealType ReduceTolerance;
};

// Inline call operators
template<bool Complex>
inline ComplexType GreensFunctionPart<Complex>::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing*RealType(2*MatsubaraNumber+1)); }

template<bool Complex>
inline ComplexType GreensFunctionPart<Complex>::operator()(ComplexType z) const {
    return Terms(z);
}

template<bool Complex>
inline ComplexType GreensFunctionPart<Complex>::of_tau(RealType tau) const {
    return Terms(tau, beta);
}

extern template class GreensFunctionPart<false>;
extern template class GreensFunctionPart<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GREENSFUNCTIONPART_H
