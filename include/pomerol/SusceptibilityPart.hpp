//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/SusceptibilityPart.hpp
/// \brief Part of a dynamical susceptibility in the Matsubara representation.
/// \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_SUSCEPTIBILITYPART_HPP
#define POMEROL_INCLUDE_SUSCEPTIBILITYPART_HPP

#include "DensityMatrixPart.hpp"
#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "MonomialOperatorPart.hpp"
#include "StatesClassification.hpp"
#include "TermList.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <complex>
#include <cstddef>
#include <limits>
#include <ostream>

namespace Pomerol {

/// \addtogroup Susc
/// @{

/// \brief Part of a dynamical susceptibility.
///
/// It includes contributions from all matrix elements of the following form,
/// \f[
///  \langle {\rm outer} | \hat A | {\rm inner} \rangle\langle {\rm inner}  | \hat B | {\rm outer} \rangle
/// \f]
/// with (inner, outer) being a certain pair of Hamiltonian's invariant subspaces.
/// The contributions are stored as terms of the Lehmann representation, i.e. as
/// fractions \f$\frac{R}{z - P}\f$ with real poles \f$P\f$ and complex residues \f$R\f$.
class SusceptibilityPart : public Thermal {
    friend class Susceptibility;

    /// Diagonal block of the Hamiltonian corresponding to the 'inner' subspace.
    HamiltonianPart const& HpartInner;
    /// Diagonal block of the Hamiltonian corresponding to the 'outer' subspace.
    HamiltonianPart const& HpartOuter;
    /// Diagonal block of the many-body density matrix corresponding to the 'inner' subspace.
    DensityMatrixPart const& DMpartInner;
    /// Diagonal block of the many-body density matrix corresponding to the 'outer' subspace.
    DensityMatrixPart const& DMpartOuter;

    /// Block of the monomial operator \f$\langle{\rm outer}|\hat A|{\rm inner}\rangle\f$.
    MonomialOperatorPart const& A;
    /// Block of the monomial operator \f$\langle{\rm outer}|\hat B|{\rm inner}\rangle\f$.
    MonomialOperatorPart const& B;

    /// A contribution to the Lehmann representation of a dynamical susceptibility,
    /// a fraction of the form \f$\frac{R}{z - P}\f$.
    struct Term {
        /// Residue at the pole (\f$ R \f$).
        ComplexType Residue = 0;
        /// Position of the pole (\f$ P \f$).
        RealType Pole = 0;

        /// Hasher for terms.
        struct Hash {
            /// Poles located within this energy spacing from each other produce the same hash value.
            double EnergySpacing;
            /// Constructor.
            /// \param[in] EnergySpacing Energy spacing.
            explicit Hash(double EnergySpacing) : EnergySpacing(EnergySpacing) {}
            /// Compute hash of a term.
            /// \param[in] t Term to compute hash for.
            std::size_t operator()(Term const& t) const { return hash_binned_real(t.Pole, EnergySpacing); }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&EnergySpacing, 1, MPI_DOUBLE, root, comm); }
        };

        /// Similarity predicate for terms.
        struct KeyEqual {
            /// Tolerance level used to compare positions of the pole.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to compare positions of the pole.
            explicit KeyEqual(double Tolerance) : Tolerance(Tolerance) {}
            /// Are terms similar?
            /// \param[in] t1 First term.
            /// \param[in] t2 Second term.
            bool operator()(Term const& t1, Term const& t2) const { return std::abs(t2.Pole - t1.Pole) <= Tolerance; }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /// Predicate: Does a term have a negligible residue?
        struct IsNegligible {
        private:
            /// Tolerance level used to detect negligible residues.
            double Tolerance;

        public:
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to detect negligible residues.
            explicit IsNegligible(double Tolerance) : Tolerance(Tolerance) {}
            /// Is term negligible?
            /// \param[in] t Term.
            /// \param[in] ToleranceDivisor Divide tolerance by this value.
            bool operator()(Term const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Residue) <= Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /// Constructor.
        /// \param[in] Residue Value of the residue \f$ R \f$.
        /// \param[in] Pole Position of the pole \f$ P \f$.
        Term(ComplexType Residue, RealType Pole);

        /// Substitute a complex frequency \f$ z \f$ into this term.
        /// \param[in] z Value of the frequency \f$ z \f$.
        ComplexType operator()(ComplexType z) const;

        /// Return the contribution to the imaginary-time Green's function made by this term.
        /// \param[in] tau Imaginary time point.
        /// \param[in] beta Inverse temperature.
        ComplexType operator()(RealType tau, RealType beta) const;

        /// In-place addition of terms (similarity of the terms is not checked).
        /// \param[in] AnotherTerm Term to add.
        Term& operator+=(Term const& AnotherTerm);
    };

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] T Term to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, Term const& T) {
        return os << T.Residue << "/(z - " << T.Pole << ")";
    }

    /// List of all terms contributing to this part.
    TermList<Term> Terms;

    /// Weight of the zero-energy pole.
    ComplexType ZeroPoleWeight = 0;

    /// Lehmann representation: Maximal distance between energy poles to be consider coinciding.
    RealType PoleResolution;
    /// Lehmann representation: Maximal magnitude of a term coefficient to be considered negligible.
    RealType CoefficientTolerance;

public:
    /// Constructor.
    /// \param[in] A Part of the monomial operator \f$\hat A\f$.
    /// \param[in] B Part of the monomial operator \f$\hat B\f$.
    /// \param[in] HpartInner Part of the Hamiltonian corresponding to the 'inner' subspace.
    /// \param[in] HpartOuter Part of the Hamiltonian corresponding to the 'outer' subspace.
    /// \param[in] DMpartInner Part of the many-body density matrix \f$\hat\rho\f$
    ///                        corresponding to the 'inner' subspace.
    /// \param[in] DMpartOuter Part of the many-body density matrix \f$\hat\rho\f$
    ///                        corresponding to the 'outer' subspace.
    /// \param[in] PoleResolution Lehmann representation: Maximal distance between energy
    ///                           poles to be consider coinciding.
    /// \param[in] CoefficientTolerance Lehmann representation: Maximal magnitude of a term
    ///                                 coefficient to be considered negligible.
    SusceptibilityPart(MonomialOperatorPart const& A,
                       MonomialOperatorPart const& B,
                       HamiltonianPart const& HpartInner,
                       HamiltonianPart const& HpartOuter,
                       DensityMatrixPart const& DMpartInner,
                       DensityMatrixPart const& DMpartOuter,
                       RealType PoleResolution,
                       RealType CoefficientTolerance);

    /// Compute the terms contributing to this part.
    void compute();

    /// Substitute a complex frequency \f$z\f$ into this part.
    /// \param[in] z Value of the frequency \f$z\f$.
    ComplexType operator()(ComplexType z) const;

    /// Substitute a bosonic Matsubara frequency \f$\omega_n\f$ into this part.
    /// \param[in] MatsubaraNumber Index of the Matsubara frequency \f$n\f$
    /// (\f$ \omega_n = 2\pi n/\beta \f$).
    ComplexType operator()(long MatsubaraNumber) const;

    /// Return the contribution to the imaginary-time susceptibility made by this part.
    /// \param[in] tau Imaginary time point.
    ComplexType of_tau(RealType tau) const;

private:
    // Implementation detail of compute().
    template <bool AComplex, bool BComplex> void computeImpl();
};

///@}

inline ComplexType SusceptibilityPart::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber));
}

inline ComplexType SusceptibilityPart::operator()(ComplexType z) const {
    // BOSON: add contribution of zero-energy pole
    ComplexType ZeroPole = std::abs(z) <= std::numeric_limits<RealType>::epsilon() ? ZeroPoleWeight * beta : 0;
    return Terms(z) + ZeroPole;
}

inline ComplexType SusceptibilityPart::of_tau(RealType tau) const {
    return Terms(tau, beta) + ZeroPoleWeight;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_SUSCEPTIBILITYPART_HPP
