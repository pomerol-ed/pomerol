//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/TwoParticleGFPart.hpp
/// \brief Part of a fermionic two-particle Matsubara Green's function.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_TWOPARTICLEGFPART_HPP
#define POMEROL_INCLUDE_TWOPARTICLEGFPART_HPP

#include "ComputableObject.hpp"
#include "DensityMatrixPart.hpp"
#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "MonomialOperatorPart.hpp"
#include "StatesClassification.hpp"
#include "TermList.hpp"
#include "Thermal.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <boost/functional/hash.hpp>

#include <array>
#include <complex>
#include <cstddef>

namespace Pomerol {

/// \addtogroup 2PGF
///@{

/// \brief Part of a fermionic two-particle Matsubara Green's function.
///
/// It includes contributions from all matrix elements of the following form,
/// \f[
///  \langle {\rm S_1}| \hat O_1 |{\rm S_2}\rangle
///  \langle {\rm S_2}| \hat O_2 |{\rm S_3} \rangle
///  \langle {\rm S_3}| \hat O_3 |{\rm S_4} \rangle
///  \langle {\rm S_4}| c^\dagger_l |{\rm S_1} \rangle,
/// \f]
/// where \f$\{\hat O_1, \hat O_2, \hat O_3\}\f$ is a permutation of operators \f$\{c_i, c_j, c^\dagger_k\}\f$
/// and \f${\rm S_1}, {\rm S_2}, {\rm S_3}, {\rm S_4}\f$ are invariant subspaces of the Hamiltonian.
/// The contributions are stored as terms of the Lehmann representation. There are two kinds of terms contributing
/// to the expansion, so called resonant (\ref ResonantTerm) and non-resonant (\ref NonResonantTerm) terms.
class TwoParticleGFPart : public Thermal, ComputableObject {

    friend class TwoParticleGF;
    friend class TwoParticleGFContainer;

public:
    /// \brief A non-resonant term in the Lehmann representation of \ref TwoParticleGF.
    ///
    /// A non-resonant term of the Lehmann representation. It is parametrized by
    /// a complex coefficient \f$C\f$ and positions of real poles \f$P_1, P_2, P_3\f$.
    /// Depending on the value of the \ref isz4 flag, an explicit expression for the term reads
    ///
    /// \li \f$\frac{C}{(z_1-P_1)(z_2-P_2)(z_3-P_3)}\f$ for \ref isz4 == false,
    /// \li \f$\frac{C}{(z_1-P_1)(z_1+z_2+z_3-P_1-P_2-P_3)(z_3-P_3)}\f$ for \ref isz4 == true.
    struct NonResonantTerm {
        /// Coefficient \f$C\f$.
        ComplexType Coeff = 0;

        /// Poles \f$P_1\f$, \f$P_2\f$, \f$P_3\f$.
        std::array<RealType, 3> Poles = {{0, 0, 0}};

        /// Are we using \f$z_4=z_1+z_2+z_3\f$ instead of \f$z_2\f$ in this term?
        bool isz4 = false;

        /// Weight \f$W\f$ used in addition of terms with different poles.
        /// \see \ref operator+=()
        long Weight = 0;

        /// Hasher for non-resonant terms.
        struct Hash {
            /// Poles located within this energy spacing from each other produce the same hash value.
            double EnergySpacing;
            /// Constructor.
            /// \param[in] EnergySpacing Energy spacing.
            explicit Hash(double EnergySpacing = 1e-8) : EnergySpacing(EnergySpacing) {}
            /// Compute hash of a term.
            /// \param[in] t Term to compute hash for.
            std::size_t operator()(NonResonantTerm const& t) const {
                auto h = boost::hash<std::tuple<bool, std::size_t, std::size_t, std::size_t>>();
                return h(std::make_tuple(t.isz4,
                                         hash_binned_real(t.Poles[0], EnergySpacing),
                                         hash_binned_real(t.Poles[1], EnergySpacing),
                                         hash_binned_real(t.Poles[2], EnergySpacing)));
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&EnergySpacing, 1, MPI_DOUBLE, root, comm); }
        };

        /// Similarity predicate for non-resonant terms.
        struct KeyEqual {
            /// Tolerance level used to compare positions of the poles.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to compare positions of the pole.
            explicit KeyEqual(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            /// Are terms similar?
            /// \param[in] t1 First term.
            /// \param[in] t2 Second term.
            bool operator()(NonResonantTerm const& t1, NonResonantTerm const& t2) const {
                return t2.isz4 == t1.isz4 && (std::abs(t2.Poles[0] - t1.Poles[0]) < Tolerance) &&
                       (std::abs(t2.Poles[1] - t1.Poles[1]) < Tolerance) &&
                       (std::abs(t2.Poles[2] - t1.Poles[2]) < Tolerance);
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /// Predicate: Does a term have a negligible residue?
        struct IsNegligible {
            /// Tolerance level used to detect negligible residues.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to detect negligible residues.
            explicit IsNegligible(double Tolerance = 1e-16) : Tolerance(Tolerance) {}
            /// Is term negligible?
            /// \param[in] t Term.
            /// \param[in] ToleranceDivisor Divide tolerance by this value.
            bool operator()(NonResonantTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Coeff) < Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        NonResonantTerm() = default;

        /// Constructor.
        /// \param[in] Coeff Coefficient of the term \f$C\f$.
        /// \param[in] P1 Pole \f$P_1\f$.
        /// \param[in] P2 Pole \f$P_2\f$.
        /// \param[in] P3 Pole \f$P_3\f$.
        /// \param[in] isz4 Are we using \f$z_4=z_1+z_2+z_3\f$ instead of \f$z_2\f$ in this term?
        inline NonResonantTerm(ComplexType Coeff, RealType P1, RealType P2, RealType P3, bool isz4)
            : Coeff(Coeff), Poles{P1, P2, P3}, isz4(isz4), Weight(1) {}

        /// Substitute complex frequencies \f$z_1, z_2, z_3\f$ into this term.
        /// \param[in] z1 Complex frequency \f$z_1\f$.
        /// \param[in] z2 Complex frequency \f$z_2\f$.
        /// \param[in] z3 Complex frequency \f$z_3\f$.
        ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;

        /// Add a non-resonant term to this term.
        ///
        /// This operator does not check similarity of the terms!
        /// Parameters of this term are updated as follows.
        /// \li Coeff += AnotherTerm.Coeff
        /// \li Poles[i] = (Poles[i] * Weight + AnotherTerm.Poles[i] * AnotherTerm.Weight) /
        ///                (Weight + AnotherTerm.Weight)
        /// \li Weight += AnotherTerm.Weight
        /// \param[in] AnotherTerm Term to add.
        NonResonantTerm& operator+=(NonResonantTerm const& AnotherTerm);

        /// Create and commit an MPI datatype for \ref NonResonantTerm.
        static MPI_Datatype mpi_datatype();
    };

    /// \brief A resonant term in the Lehmann representation of \ref TwoParticleGF.
    ///
    /// It is parametrized by
    /// two complex coefficients \f$R\f$ and \f$N\f$, and positions of real poles \f$P_1, P_2, P_3\f$.
    /// Depending on the value of the \ref isz1z2 flag, an explicit expression for the term reads
    /// \li \f$
    ///   \frac{1}{(z_1-P_1)(z_3-P_3)}
    ///   \left( R \delta(z_1+z_2-P_1-P_2) + N \frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2} \right)
    /// \f$ for \ref isz1z2 == true,
    /// \li \f$
    ///   \frac{1}{(z_1-P_1)(z_3-P_3)}
    ///   \left( R \delta(z_2+z_3-P_2-P_3) + N \frac{1 - \delta(z_2+z_3-P_2-P_3)}{z_2+z_3-P_2-P_3} \right)
    /// \f$ for \ref isz1z2 == false.
    struct ResonantTerm {

        /// Coefficient \f$R\f$.
        ComplexType ResCoeff = 0;
        /// Coefficient \f$N\f$.
        ComplexType NonResCoeff = 0;

        /// Poles \f$P_1\f$, \f$P_2\f$, \f$P_3\f$.
        std::array<RealType, 3> Poles = {{0, 0, 0}};

        /// Are we using \f$ \delta(z_1+z_2-P_1-P_2) \f$ resonance condition?
        /// If not, we are using \f$ \delta(z_2+z_3-P_2-P_3) \f$.
        bool isz1z2 = false;

        /// Weight \f$W\f$ used in addition of terms with different poles.
        /// \see \ref operator+=()
        long Weight = 0;

        /// Hasher for resonant terms.
        struct Hash {
            /// Poles located within this energy spacing from each other produce the same hash value.
            double EnergySpacing;
            /// Constructor.
            /// \param[in] EnergySpacing Energy spacing.
            explicit Hash(double EnergySpacing = 1e-8) : EnergySpacing(EnergySpacing) {}
            /// Compute hash of a term.
            /// \param[in] t Term to compute hash for.
            std::size_t operator()(ResonantTerm const& t) const {
                auto h = boost::hash<std::tuple<bool, std::size_t, std::size_t, std::size_t>>();
                return h(std::make_tuple(t.isz1z2,
                                         hash_binned_real(t.Poles[0], EnergySpacing),
                                         hash_binned_real(t.Poles[1], EnergySpacing),
                                         hash_binned_real(t.Poles[2], EnergySpacing)));
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&EnergySpacing, 1, MPI_DOUBLE, root, comm); }
        };

        /// Similarity predicate for resonant terms.
        struct KeyEqual {
            /// Tolerance level used to compare positions of the poles.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to compare positions of the pole.
            explicit KeyEqual(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            /// Are terms similar?
            /// \param[in] t1 First term.
            /// \param[in] t2 Second term.
            bool operator()(ResonantTerm const& t1, ResonantTerm const& t2) const {
                return t2.isz1z2 == t1.isz1z2 && (std::abs(t2.Poles[0] - t1.Poles[0]) < Tolerance) &&
                       (std::abs(t2.Poles[1] - t1.Poles[1]) < Tolerance) &&
                       (std::abs(t2.Poles[2] - t1.Poles[2]) < Tolerance);
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        /// Predicate: Does a term have a negligible residue?
        struct IsNegligible {
            /// Tolerance level used to detect negligible residues.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to detect negligible residues.
            explicit IsNegligible(double Tolerance = 1e-16) : Tolerance(Tolerance) {}
            /// Is term negligible?
            /// \param[in] t Term.
            /// \param[in] ToleranceDivisor Divide tolerance by this value.
            bool operator()(ResonantTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.ResCoeff) < Tolerance / ToleranceDivisor &&
                       std::abs(t.NonResCoeff) < Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        ResonantTerm() = default;

        /// Constructor.
        /// \param[in] ResCoeff Numerator of the term for the resonant case, \f$R\f$.
        /// \param[in] NonResCoeff Numerator of the term for the non-resonant case, \f$N\f$.
        /// \param[in] P1 Pole \f$P_1\f$.
        /// \param[in] P2 Pole \f$P_2\f$.
        /// \param[in] P3 Pole \f$P_3\f$.
        /// \param[in] isz1z2 Are we using the \f$\delta(z_1+z_2-P_1-P_2)\f$ resonance condition?
        inline ResonantTerm(ComplexType ResCoeff,
                            ComplexType NonResCoeff,
                            RealType P1,
                            RealType P2,
                            RealType P3,
                            bool isz1z2)
            : ResCoeff(ResCoeff), NonResCoeff(NonResCoeff), Poles{P1, P2, P3}, isz1z2(isz1z2), Weight(1) {}

        /// Substitute complex frequencies \f$z_1, z_2, z_3\f$ into this term.
        /// \param[in] z1 Complex frequency \f$z_1\f$.
        /// \param[in] z2 Complex frequency \f$z_2\f$.
        /// \param[in] z3 Complex frequency \f$z_3\f$.
        /// \param[in] DeltaTolerance Tolerance for the resonance detection.
        ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3, RealType DeltaTolerance = 1e-16) const;

        /// Add a resonant term to this term.
        ///
        /// This operator does not check similarity of the terms!
        /// Parameters of this term are updated as follows.
        /// \li ResCoeff += AnotherTerm.ResCoeff
        /// \li NonResCoeff += AnotherTerm.NonResCoeff
        /// \li Poles[i] = (Poles[i] * Weight + AnotherTerm.Poles[i] * AnotherTerm.Weight) /
        ///                (Weight + AnotherTerm.Weight)
        /// \li Weight += AnotherTerm.Weight
        /// \param[in] AnotherTerm Term to add.
        ResonantTerm& operator+=(ResonantTerm const& AnotherTerm);

        /// Create and commit an MPI datatype for \ref ResonantTerm.
        static MPI_Datatype mpi_datatype();
    };

private:
    /// Part of the field operator \f$\hat O_1\f$.
    MonomialOperatorPart const& O1;
    /// Part of the field operator \f$\hat O_2\f$.
    MonomialOperatorPart const& O2;
    /// Part of the field operator \f$\hat O_3\f$.
    MonomialOperatorPart const& O3;
    /// Part of the creation operator \f$\hat c^\dagger_l\f$.
    MonomialOperatorPart const& CX4;

    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_1}\f$.
    HamiltonianPart const& Hpart1;
    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_2}\f$.
    HamiltonianPart const& Hpart2;
    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_3}\f$.
    HamiltonianPart const& Hpart3;
    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_4}\f$.
    HamiltonianPart const& Hpart4;

    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_1}\f$.
    DensityMatrixPart const& DMpart1;
    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_2}\f$.
    DensityMatrixPart const& DMpart2;
    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_3}\f$.
    DensityMatrixPart const& DMpart3;
    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_4}\f$.
    DensityMatrixPart const& DMpart4;

    /// Permutation of the operators \f$\{c_i, c_j, c^\dagger_k\}\f$ for this part.
    Permutation3 Permutation;

    /// List of all non-resonant terms contributing to this part.
    TermList<NonResonantTerm> NonResonantTerms;
    /// List of all resonant terms contributing to this part.
    TermList<ResonantTerm> ResonantTerms;

    /// Adds a multi-term that has the following form:
    /// \f[
    /// \frac{1}{(z_1-P_1)(z_3-P_3)}
    ///         \left(\frac{C_4}{z_1+z_2+z_3-P_1-P_2-P_3} + \frac{C_2}{z_2-P_2} \right. +
    /// \f]
    /// \f[     \left.
    ///         + R_{12}\delta(z_1+z_2-P_1-P_2)
    ///         + N_{12}\frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2}
    ///         + R_{23}\delta(z_2+z_3-P_2-P_3)
    ///         + N_{23}\frac{1 - \delta(z_2+z_3-P_2-P_3)}{z_2+z_3-P_2-P_3}
    ///         \right),
    /// \f]
    /// where
    /// \f{eqnarray*}{
    ///      P_1 = E_j - E_i \\
    ///      P_2 = E_k - E_j \\
    ///      P_3 = E_l - E_k \\
    ///      C_2 = -C(w_j + w_k) \\
    ///      C_4 = C(w_i + w_l) \\
    ///      R_{12} = C\beta w_i \\
    ///      N_{12} = C(w_k - w_i) \\
    ///      R_{23} = -C\beta w_j \\
    ///      N_{23} = C(w_j - w_l)
    /// \f}
    ///
    /// In fact this is a slightly rewritten form of an equation for \f$\phi\f$ from
    /// <em>H. Hafermann et al 2009 EPL 85 27007</em>.
    ///
    /// \param[in] Coeff Common prefactor \f$C\f$ for coefficients \f$C_2\f$, \f$C_4\f$,
    ///              \f$R_{12}\f$, \f$N_{12}\f$, \f$R_{23}\f$, \f$N_{23}\f$.
    /// \param[in] beta Inverse temperature.
    /// \param[in] Ei The first energy level \f$E_i\f$.
    /// \param[in] Ej The second energy level \f$E_j\f$.
    /// \param[in] Ek The third energy level \f$E_k\f$.
    /// \param[in] El The fourth energy level \f$E_l\f$.
    /// \param[in] Wi The first weight \f$w_i\f$.
    /// \param[in] Wj The second weight \f$w_j\f$.
    /// \param[in] Wk The third weight \f$w_k\f$.
    /// \param[in] Wl The fourth weight \f$w_l\f$.
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

    /// A difference in energies with magnitude below this value is treated as zero.
    RealType ReduceResonanceTolerance = 1e-8;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account.
    RealType CoefficientTolerance = 1e-16;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account with respect to
    /// the amount of terms.
    RealType MultiTermCoefficientTolerance = 1e-5;

    // compute() implementation details.
    template <bool Complex> void computeImpl();

public:
    /// Constructor.
    /// \param[in] O1 Part of the field operator \f$\hat O_1\f$.
    /// \param[in] O2 Part of the field operator \f$\hat O_2\f$.
    /// \param[in] O3 Part of the field operator \f$\hat O_3\f$.
    /// \param[in] CX4 Part of the creation operator \f$\hat c^\dagger_l\f$.
    /// \param[in] Hpart1 Part of the Hamiltonian corresponding to the subspace \f${\rm S_1}\f$.
    /// \param[in] Hpart2 Part of the Hamiltonian corresponding to the subspace \f${\rm S_2}\f$.
    /// \param[in] Hpart3 Part of the Hamiltonian corresponding to the subspace \f${\rm S_3}\f$.
    /// \param[in] Hpart4 Part of the Hamiltonian corresponding to the subspace \f${\rm S_4}\f$.
    /// \param[in] DMpart1 Part of the many-body density matrix corresponding to the subspace \f${\rm S_1}\f$.
    /// \param[in] DMpart2 Part of the many-body density matrix corresponding to the subspace \f${\rm S_2}\f$.
    /// \param[in] DMpart3 Part of the many-body density matrix corresponding to the subspace \f${\rm S_3}\f$.
    /// \param[in] DMpart4 Part of the many-body density matrix corresponding to the subspace \f${\rm S_4}\f$.
    /// \param[in] Permutation Permutation of operators \f$\{c_i, c_j, c^\dagger_k\}\f$ for this part.
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

    /// Compute the terms contributing to this part.
    void compute();

    /// Purge all terms.
    void clear();

    /// Substitute complex frequencies \f$z_1, z_2, z_3\f$ into this part.
    /// \param[in] z1 First frequency \f$z_1\f$.
    /// \param[in] z2 Second frequency \f$z_2\f$.
    /// \param[in] z3 Third frequency \f$z_3\f$.
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;
    /// Substitute Matsubara frequencies \f$i\omega_{n_1}, i\omega_{n_2}, i\omega_{n_3}\f$ into this part.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    /// \param[in] MatsubaraNumber3 Index of the third Matsubara frequency
    ///                             \f$n_3\f$ (\f$\omega_{n_3}=\pi(2n_3+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /// Return the number of resonant terms.
    std::size_t getNumResonantTerms() const { return ResonantTerms.size(); }
    /// Return the number of non-resonant terms.
    std::size_t getNumNonResonantTerms() const { return NonResonantTerms.size(); }

    /// Return the permutation of operators \f$\{c_i, c_j, c^\dagger_k\}\f$ for this part.
    Permutation3 const& getPermutation() const { return Permutation; }

    /// Access the list of the resonant terms.
    TermList<TwoParticleGFPart::ResonantTerm> const& getResonantTerms() const { return ResonantTerms; }
    /// Access the list of the non-resonant terms.
    TermList<TwoParticleGFPart::NonResonantTerm> const& getNonResonantTerms() const { return NonResonantTerms; }
};

///@}

inline ComplexType
TwoParticleGFPart::NonResonantTerm::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const {
    return isz4 ? Coeff / ((z1 - Poles[0]) * (z1 + z2 + z3 - Poles[0] - Poles[1] - Poles[2]) * (z3 - Poles[2])) :
                  Coeff / ((z1 - Poles[0]) * (z2 - Poles[1]) * (z3 - Poles[2]));
}

inline ComplexType TwoParticleGFPart::ResonantTerm::operator()(ComplexType z1,
                                                               ComplexType z2,
                                                               ComplexType z3,
                                                               RealType DeltaTolerance) const {
    ComplexType Diff;
    if(isz1z2) {
        Diff = z1 + z2 - Poles[0] - Poles[1];
        return (std::abs(Diff) < DeltaTolerance ? ResCoeff : (NonResCoeff / Diff)) /
               ((z1 - Poles[0]) * (z3 - Poles[2]));
    } else {
        Diff = z2 + z3 - Poles[1] - Poles[2];
        return (std::abs(Diff) < DeltaTolerance ? ResCoeff : (NonResCoeff / Diff)) /
               ((z1 - Poles[0]) * (z3 - Poles[2]));
    }
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGFPART_HPP
