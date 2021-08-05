/** \file include/pomerol/GreensFunctionPart.h
** \brief Part of a Green's function for a given set of quantum numbers.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_GREENSFUNCTIONPART_H
#define POMEROL_INCLUDE_POMEROL_GREENSFUNCTIONPART_H

#include "DensityMatrixPart.hpp"
#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "MonomialOperatorPart.hpp"
#include "StatesClassification.hpp"
#include "TermList.hpp"
#include "Thermal.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <complex>
#include <cstddef>
#include <ostream>

namespace Pomerol {

/** This class represents a part of a Green's function.
 * Every part describes all transitions allowed by selection rules
 * between a given pair of Hamiltonian blocks.
 */
class GreensFunctionPart : public Thermal {
    /** A reference to a part of a Hamiltonian (inner index iterates through it). */
    HamiltonianPart const& HpartInner;
    /** A reference to a part of a Hamiltonian (outer index iterates through it). */
    HamiltonianPart const& HpartOuter;
    /** A reference to a part of a density matrix (the part corresponding to HpartInner). */
    DensityMatrixPart const& DMpartInner;
    /** A reference to a part of a density matrix (the part corresponding to HpartOuter). */
    DensityMatrixPart const& DMpartOuter;

    /** A reference to a part of an annihilation operator. */
    MonomialOperatorPart const& C;
    /** A reference to a part of a creation operator. */
    MonomialOperatorPart const& CX;

    /** Every term is a fraction \f$ \frac{R}{z - P} \f$. */
    struct Term {
        /** Residue at the pole (\f$ R \f$). */
        ComplexType Residue;
        /** Position of the pole (\f$ P \f$). */
        RealType Pole;

        /** Comparator object for terms */
        struct Compare {
            double const Tolerance;
            Compare(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            bool operator()(Term const& t1, Term const& t2) const { return t2.Pole - t1.Pole >= Tolerance; }
        };

        /** Does term have a negligible residue? */
        struct IsNegligible {
            double Tolerance;
            IsNegligible(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            bool operator()(Term const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Residue) < Tolerance / ToleranceDivisor;
            }
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
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
        Term& operator+=(Term const& AnotherTerm);
    };
    /** A stream insertion operator for type GreensTerm.
     * \param[in] out An output stream to insert to.
     * \param[in] Term A term to be inserted.
     */
    friend std::ostream& operator<<(std::ostream& os, Term const& T) {
        return os << T.Residue << "/(z - " << T.Pole << ")";
    }

    /** A list of all terms. */
    TermList<Term> Terms;

    /** A matrix element with magnitude less than this value is treated as zero. */
    RealType const MatrixElementTolerance = 1e-8;

public:
    /** Constructor.
     * \param[in] C A reference to a part of an annihilation operator.
     * \param[in] CX A reference to a part of a creation operator.
     * \param[in] HpartInner A reference to a part of the Hamiltonian (inner index).
     * \param[in] HpartOuter A reference to a part of the Hamiltonian (outer index).
     * \param[in] DMpartInner A reference to a part of the density matrix (inner index).
     * \param[in] DMpartOuter A reference to a part of the density matrix (outer index).
     */
    GreensFunctionPart(MonomialOperatorPart const& C,
                       MonomialOperatorPart const& CX,
                       HamiltonianPart const& HpartInner,
                       HamiltonianPart const& HpartOuter,
                       DensityMatrixPart const& DMpartInner,
                       DensityMatrixPart const& DMpartOuter);

    /** Iterates over all matrix elements and fills the list of terms. */
    void compute();

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
    RealType const ReduceResonanceTolerance = 1e-8;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
    RealType const ReduceTolerance = 1e-8;

private:
    template <bool Complex> void computeImpl();
};

// Inline call operators
inline ComplexType GreensFunctionPart::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber + 1));
}

inline ComplexType GreensFunctionPart::operator()(ComplexType z) const {
    return Terms(z);
}

inline ComplexType GreensFunctionPart::of_tau(RealType tau) const {
    return Terms(tau, beta);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_GREENSFUNCTIONPART_H
