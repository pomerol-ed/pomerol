//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/TwoParticleGFPart.h
** \brief Part of a two-particle Green's function.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_TWOPARTICLEGFPART_H
#define POMEROL_INCLUDE_TWOPARTICLEGFPART_H

#include "ComputableObject.hpp"
#include "DensityMatrixPart.hpp"
#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "MonomialOperatorPart.hpp"
#include "StatesClassification.hpp"
#include "TermList.hpp"
#include "Thermal.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <array>
#include <complex>
#include <cstddef>

namespace Pomerol {

/** This class represents a part of a two-particle Green's function.
 * Every part describes one 'world stripe' of four operators.
 */
class TwoParticleGFPart : public Thermal, ComputableObject {

    friend class TwoParticleGF;
    friend class TwoParticleGFContainer;

public:
    /** A non-resonant term has the following form:
     * \f[
     * \frac{C}{(z_1-P_1)(z_2-P_2)(z_3-P_3)}
     * \f]
     * if isz4 == false, and
     * \f[
     * \frac{C}{(z_1-P_1)(z_1+z_2+z_3-P_1-P_2-P_3)(z_3-P_3)}
     * \f]
     * otherwise.
     */
    struct NonResonantTerm {
        /** Coefficient \f$ C \f$. */
        ComplexType Coeff;

        /** Poles \f$ P_1 \f$, \f$ P_2 \f$, \f$ P_3 \f$. */
        std::array<RealType, 3> Poles;

        /** Are we using \f$ z_4 \f$ instead of \f$ z_2 \f$ in this term? */
        bool isz4;

        /** A statistical weight of current term for averaging ( when averaging formula (*this.weight + other.weight)/(*this.weight+other.weight) is used */
        long Weight;

        /** Comparator object for terms */
        struct Compare {
            double Tolerance;
            Compare(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            bool real_eq(RealType x1, RealType x2) const { return std::abs(x1 - x2) < Tolerance; }
            bool operator()(NonResonantTerm const& t1, NonResonantTerm const& t2) const {
                if(t1.isz4 == t2.isz4) {
                    return !real_eq(t1.Poles[0], t2.Poles[0]) ?
                               t1.Poles[0] < t2.Poles[0] :
                               (!real_eq(t1.Poles[1], t2.Poles[1]) ? t1.Poles[1] < t2.Poles[1] :
                                                                     (t2.Poles[2] - t1.Poles[2] >= Tolerance));
                } else
                    return t1.isz4 < t2.isz4;
            }
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /** Does term have a negligible residue? */
        struct IsNegligible {
            double Tolerance;
            IsNegligible(double Tolerance = 1e-16) : Tolerance(Tolerance) {}
            bool operator()(NonResonantTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Coeff) < Tolerance / ToleranceDivisor;
            }
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        NonResonantTerm() = default;

        /** Constructor.
        * \param[in] Coeff Numerator of the term.
        * \param[in] P1 Pole P1.
        * \param[in] P2 Pole P2.
        * \param[in] P3 Pole P3.
        * \param[in] isz4 Are we using \f$ z_4 \f$ instead of \f$ z_2 \f$ in this term?
        */
        inline NonResonantTerm(ComplexType Coeff, RealType P1, RealType P2, RealType P3, bool isz4)
            : Coeff(Coeff), isz4(isz4) {
            Poles[0] = P1;
            Poles[1] = P2;
            Poles[2] = P3;
            Weight = 1;
        }

        /** Returns a contribution to the two-particle Green's function made by this term.
        * \param[in] z1 Complex frequency \f$ z_1 \f$.
        * \param[in] z2 Complex frequency \f$ z_2 \f$.
        * \param[in] z3 Complex frequency \f$ z_3 \f$.
        */
        ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;

        /** This operator add a non-resonant term to this one.
        * It does not check the similarity of the terms!
        * \param[in] AnotherTerm Another term to add to this.
        */
        NonResonantTerm& operator+=(NonResonantTerm const& AnotherTerm);

        /** Create and commit an MPI datatype for NonResonantTerm */
        static MPI_Datatype mpi_datatype();
    };

    /** A resonant term has the following form:
     * \f[
     * \frac{1}{(z_1-P_1)(z_3-P_3)}
     *   \left( R \delta(z_1+z_2-P_1-P_2) + N \frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2} \right)
     * \f]
     */

    struct ResonantTerm {

        /** Coefficient \f$ R \f$. */
        ComplexType ResCoeff;
        /** Coefficient \f$ N \f$. */
        ComplexType NonResCoeff;

        /** Poles \f$ P_1 \f$, \f$ P_2 \f$, \f$ P_3 \f$. */
        std::array<RealType, 3> Poles;

        /** Are we using \f$ \delta(z_1+z_2-P_1-P_2) \f$ resonance condition?
        Otherwise we are using \f$ \delta(z_2+z_3-P_2-P_3) \f$. */
        bool isz1z2;

        /** A statistical weight of current term for averaging ( when averaging formula (*this.weight + other.weight)/(*this.weight+other.weight) is used */
        long Weight;

        /** Comparator object for terms */
        struct Compare {
            double Tolerance;
            Compare(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            bool real_eq(RealType x1, RealType x2) const { return std::abs(x1 - x2) < Tolerance; }
            bool operator()(ResonantTerm const& t1, ResonantTerm const& t2) const {
                if(t1.isz1z2 == t2.isz1z2) {
                    return !real_eq(t1.Poles[0], t2.Poles[0]) ?
                               t1.Poles[0] < t2.Poles[0] :
                               (!real_eq(t1.Poles[1], t2.Poles[1]) ? t1.Poles[1] < t2.Poles[1] :
                                                                     (t2.Poles[2] - t1.Poles[2] >= Tolerance));
                } else
                    return t1.isz1z2 < t2.isz1z2;
            }
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /** Does term have a negligible residue? */
        struct IsNegligible {
            double Tolerance;
            IsNegligible(double Tolerance = 1e-16) : Tolerance(Tolerance) {}
            bool operator()(ResonantTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.ResCoeff) < Tolerance / ToleranceDivisor &&
                       std::abs(t.NonResCoeff) < Tolerance / ToleranceDivisor;
            }
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        ResonantTerm() = default;

        /** Constructor.
        * \param[in] ResCoeff Numerator of the term for a resonant case.
        * \param[in] NonResCoeff Numerator of the term for a non-resonant case.
        * \param[in] P1 Pole P1.
        * \param[in] P2 Pole P2.
        * \param[in] P3 Pole P3.
        * \param[in] isz1z2 Are we using \f$ \delta(z_1+z_2-P_1-P_2) \f$ resonance condition?
        */
        inline ResonantTerm(ComplexType ResCoeff,
                            ComplexType NonResCoeff,
                            RealType P1,
                            RealType P2,
                            RealType P3,
                            bool isz1z2)
            : ResCoeff(ResCoeff), NonResCoeff(NonResCoeff), isz1z2(isz1z2) {
            Poles[0] = P1;
            Poles[1] = P2;
            Poles[2] = P3;
            Weight = 1;
        }

        /** Returns a contribution to the two-particle Green's function made by this term.
        * \param[in] z1 Complex frequency \f$ z_1 \f$.
        * \param[in] z2 Complex frequency \f$ z_2 \f$.
        * \param[in] z3 Complex frequency \f$ z_3 \f$.
        */
        ComplexType
        operator()(ComplexType z1, ComplexType z2, ComplexType z3, RealType KroneckerSymbolTolerance = 1e-16) const;

        /** This operator add a non-resonant term to this one.
        * It does not check the similarity of the terms!
        * \param[in] AnotherTerm Another term to add to this.
        */
        ResonantTerm& operator+=(ResonantTerm const& AnotherTerm);

        /** Create and commit an MPI datatype for ResonantTerm */
        static MPI_Datatype mpi_datatype();
    };

private:
    /** A reference to a part of the first operator. */
    MonomialOperatorPart const& O1;
    /** A reference to a part of the second operator. */
    MonomialOperatorPart const& O2;
    /** A reference to a part of the third operator. */
    MonomialOperatorPart const& O3;
    /** A reference to a part of the fourth (creation) operator. */
    MonomialOperatorPart const& CX4;

    /** A reference to the first part of a Hamiltonian. */
    HamiltonianPart const& Hpart1;
    /** A reference to the second part of a Hamiltonian. */
    HamiltonianPart const& Hpart2;
    /** A reference to the third part of a Hamiltonian. */
    HamiltonianPart const& Hpart3;
    /** A reference to the fourth part of a Hamiltonian. */
    HamiltonianPart const& Hpart4;

    /** A reference to the first part of a density matrix (the part corresponding to Hpart1). */
    DensityMatrixPart const& DMpart1;
    /** A reference to the second part of a density matrix (the part corresponding to Hpart2). */
    DensityMatrixPart const& DMpart2;
    /** A reference to the third part of a density matrix (the part corresponding to Hpart3). */
    DensityMatrixPart const& DMpart3;
    /** A reference to the fourth part of a density matrix (the part corresponding to Hpart4). */
    DensityMatrixPart const& DMpart4;

    /** A permutation of the operators for this part. */
    Permutation3 Permutation;

    /** A list of non-resonant terms. */
    TermList<NonResonantTerm> NonResonantTerms;
    /** A list of resonant terms. */
    TermList<ResonantTerm> ResonantTerms;

    /** Adds a multi-term that has the following form:
    * \f[
    * \frac{1}{(z_1-P_1)(z_3-P_3)}
    *         \left(\frac{C_4}{z_1+z_2+z_3-P_1-P_2-P_3} + \frac{C_2}{z_2-P_2} \right. +
    * \f]
    * \f[     \left.
    *         + R_{12}\delta(z_1+z_2-P_1-P_2)
    *         + N_{12}\frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2}
    *         + R_{23}\delta(z_2+z_3-P_2-P_3)
    *         + N_{23}\frac{1 - \delta(z_2+z_3-P_2-P_3)}{z_2+z_3-P_2-P_3}
    *         \right)
    * \f]
    *
    * Where
    * \f{eqnarray*}{
    *      P_1 = E_j - E_i \\
    *      P_2 = E_k - E_j \\
    *      P_3 = E_l - E_k \\
    *      C_2 = -C(w_j + w_k) \\
    *      C_4 = C(w_i + w_l) \\
    *      R_{12} = C\beta w_i \\
    *      N_{12} = C(w_k - w_i) \\
    *      R_{23} = -C\beta w_j \\
    *      N_{23} = C(w_j - w_l)
    * \f}
    *
    * In fact this is a slightly rewritten form of an equation for \f$ \phi \f$ from
    * <em>H. Hafermann et al 2009 EPL 85 27007</em>.
    *
    * \param[in] Coeff Common prefactor \f$ C \f$ for coefficients \f$ C_2 \f$, \f$ C_4 \f$,
    *              \f$ R_{12} \f$, \f$ N_{12} \f$, \f$ R_{23} \f$, \f$ N_{23} \f$.
    * \param[in] beta The inverse temperature.
    * \param[in] Ei The first energy level \f$ E_i \f$.
    * \param[in] Ej The second energy level \f$ E_j \f$.
    * \param[in] Ek The third energy level \f$ E_k \f$.
    * \param[in] El The fourth energy level \f$ E_l \f$.
    * \param[in] Wi The first weight \f$ w_i \f$.
    * \param[in] Wj The second weight \f$ w_j \f$.
    * \param[in] Wk The third weight \f$ w_k \f$.
    * \param[in] Wl The fourth weight \f$ w_l \f$.
    * \param[in] Permutation A reference to a permutation of operators for this part.
    */
    void addMultiterm(ComplexType Coeff,
                      RealType beta,
                      RealType Ei,
                      RealType Ej,
                      RealType Ek,
                      RealType El,
                      RealType Wi,
                      RealType Wj,
                      RealType Wk,
                      RealType Wl);

    /** A difference in energies with magnitude less than this value is treated as zero. default = 1e-8. */
    RealType ReduceResonanceTolerance = 1e-8;
    /** Minimal magnitude of the coefficient of a term to take it into account. default = 1e-16. */
    RealType CoefficientTolerance = 1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. default = 1e-5. */
    RealType MultiTermCoefficientTolerance = 1e-5;

    template <bool Complex> void computeImpl();

public:
    /** Constructor.
     * \param[in] O1 A reference to a part of the first operator.
     * \param[in] O2 A reference to a part of the second operator.
     * \param[in] O3 A reference to a part of the third operator.
     * \param[in] CX4 A reference to a part of the fourth (creation) operator.
     * \param[in] Hpart1 A reference to the first part of a Hamiltonian.
     * \param[in] Hpart2 A reference to the second part of a Hamiltonian.
     * \param[in] Hpart3 A reference to the third part of a Hamiltonian.
     * \param[in] Hpart4 A reference to the fourth part of a Hamiltonian.
     * \param[in] DMpart1 A reference to the first part of a density matrix.
     * \param[in] DMpart2 A reference to the second part of a density matrix.
     * \param[in] DMpart3 A reference to the third part of a density matrix.
     * \param[in] DMpart4 A reference to the fourth part of a density matrix.
     * \param[in] Permutation A permutation of the operators for this part.
     */
    TwoParticleGFPart(MonomialOperatorPart const& O1,
                      MonomialOperatorPart const& O2,
                      MonomialOperatorPart const& O3,
                      MonomialOperatorPart const& CX4,
                      HamiltonianPart const& Hpart1,
                      HamiltonianPart const& Hpart2,
                      HamiltonianPart const& Hpart3,
                      HamiltonianPart const& Hpart4,
                      DensityMatrixPart const& DMpart1,
                      DensityMatrixPart const& DMpart2,
                      DensityMatrixPart const& DMpart3,
                      DensityMatrixPart const& DMpart4,
                      Permutation3 Permutation);

    /** Actually computes the part. */
    void compute();

    /** Purges all terms. */
    void clear();

    /** Returns the value of the Green's function calculated at a given frequency (ignores precomputed values).
    * \param[in] z1 Frequency 1
    * \param[in] z2 Frequency 2
    * \param[in] z3 Frequency 3
    */
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;
    /** Returns a contribution to the two-particle Green's function made by this part.
    * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
    * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
    * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
    */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /** Returns the number of resonant terms in the cache. */
    std::size_t getNumResonantTerms() const { return ResonantTerms.size(); }
    /** Returns the number of non-resonant terms in the cache. */
    std::size_t getNumNonResonantTerms() const { return NonResonantTerms.size(); }

    /** Returns a Permutation3 of the current part */
    Permutation3 const& getPermutation() const { return Permutation; }

    /** Return the list of Resonant Terms */
    TermList<TwoParticleGFPart::ResonantTerm> const& getResonantTerms() const { return ResonantTerms; }
    /** Return the list of NonResonantTerms */
    TermList<TwoParticleGFPart::NonResonantTerm> const& getNonResonantTerms() const { return NonResonantTerms; }
};

inline ComplexType
TwoParticleGFPart::NonResonantTerm::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const {
    return isz4 ? Coeff / ((z1 - Poles[0]) * (z1 + z2 + z3 - Poles[0] - Poles[1] - Poles[2]) * (z3 - Poles[2])) :
                  Coeff / ((z1 - Poles[0]) * (z2 - Poles[1]) * (z3 - Poles[2]));
}

inline ComplexType TwoParticleGFPart::ResonantTerm::operator()(ComplexType z1,
                                                               ComplexType z2,
                                                               ComplexType z3,
                                                               RealType KroneckerSymbolTolerance) const {
    ComplexType Diff;
    if(isz1z2) {
        Diff = z1 + z2 - Poles[0] - Poles[1];
        return (std::abs(Diff) < KroneckerSymbolTolerance ? ResCoeff : (NonResCoeff / Diff)) /
               ((z1 - Poles[0]) * (z3 - Poles[2]));
    } else {
        Diff = z2 + z3 - Poles[1] - Poles[2];
        return (std::abs(Diff) < KroneckerSymbolTolerance ? ResCoeff : (NonResCoeff / Diff)) /
               ((z1 - Poles[0]) * (z3 - Poles[2]));
    }
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGFPART_H
