//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/ThreePointSusceptibilityPart.hpp
/// \brief Part of a 3-point susceptibility in the Matsubara representation.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYPART_HPP
#define POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYPART_HPP

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
#include <tuple>

namespace Pomerol {

/// \addtogroup 3PSusc
///@{

/// \brief Part of a 3-point susceptibility.
///
/// It includes contributions from all matrix elements of the following form,
/// \f[
///  \langle {\rm S_1}| \hat F_1 |{\rm S_2}\rangle
///  \langle {\rm S_2}| \hat F_2 |{\rm S_3} \rangle
///  \langle {\rm S_3}| \hat B_1 |{\rm S_4} \rangle
///  \langle {\rm S_4}| \hat B_2 |{\rm S_1} \rangle,
/// \f]
/// where \f$\hat F_1, \hat F_2, \hat B_1, \hat B_2\f$ are fermionic field operators,
/// and \f$\hat B = \hat B_1 \hat B_2\f$.
/// \f${\rm S_1}, {\rm S_2}, {\rm S_3}, {\rm S_4}\f$ are invariant subspaces of the Hamiltonian.
/// The contributions are stored as terms of the Lehmann representation. There are three kinds of terms contributing
/// to the expansion, so called resonant (\ref ResonantTerm), non-resonant fermion-fermion (\ref NonResonantFFTerm)
/// and non-resonant fermion-boson (\ref NonResonantFBTerm) terms.
class ThreePointSusceptibilityPart : public Thermal, ComputableObject {

    friend class ThreePointSusceptibility;
    friend class ThreePointSusceptibilityContainer;

public:
    /// \brief A non-resonant fermion-fermion term in the Lehmann representation of \ref ThreePointSusceptibility.
    ///
    /// It is parametrized by a complex coefficient \f$C\f$ and positions of real poles \f$P_1, P_2\f$.
    /// An explicit expression for the term reads \f$\frac{C}{(z_1-P_1)(z_2-P_2)}\f$.
    struct NonResonantFFTerm {
        /// Coefficient \f$C\f$.
        ComplexType Coeff = 0;

        /// Poles \f$P_1\f$, \f$P_2\f$.
        std::array<RealType, 2> Poles = {{0, 0}};

        /// Weight \f$W\f$ used in addition of terms with different poles.
        /// \see \ref operator+=()
        long Weight = 0;

        /// Hasher for non-resonant fermion-fermion terms.
        struct Hash {
            /// Poles located within this energy spacing from each other produce the same hash value.
            double EnergySpacing;
            /// Constructor.
            /// \param[in] EnergySpacing Energy spacing.
            explicit Hash(double EnergySpacing = 1e-8) : EnergySpacing(EnergySpacing) {}
            /// Compute hash of a term.
            /// \param[in] t Term to compute hash for.
            std::size_t operator()(NonResonantFFTerm const& t) const {
                auto h = boost::hash<std::tuple<std::size_t, std::size_t>>();
                return h(std::make_tuple(hash_binned_real(t.Poles[0], EnergySpacing),
                                         hash_binned_real(t.Poles[1], EnergySpacing)));
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&EnergySpacing, 1, MPI_DOUBLE, root, comm); }
        };

        /// Similarity predicate for non-resonant fermion-fermion terms.
        struct KeyEqual {
            /// Tolerance level used to compare positions of the poles.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to compare positions of the pole.
            explicit KeyEqual(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            /// Are terms similar?
            /// \param[in] t1 First term.
            /// \param[in] t2 Second term.
            bool operator()(NonResonantFFTerm const& t1, NonResonantFFTerm const& t2) const {
                return (std::abs(t2.Poles[0] - t1.Poles[0]) < Tolerance) &&
                       (std::abs(t2.Poles[1] - t1.Poles[1]) < Tolerance);
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
            bool operator()(NonResonantFFTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Coeff) < Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        NonResonantFFTerm() = default;

        /// Constructor.
        /// \param[in] Coeff Coefficient of the term \f$C\f$.
        /// \param[in] P1 Pole \f$P_1\f$.
        /// \param[in] P2 Pole \f$P_2\f$.
        inline NonResonantFFTerm(ComplexType Coeff, RealType P1, RealType P2)
            : Coeff(Coeff), Poles{P1, P2}, Weight(1) {}

        /// Substitute complex frequencies \f$z_1, z_2\f$ into this term.
        /// \param[in] z1 Complex frequency \f$z_1\f$.
        /// \param[in] z2 Complex frequency \f$z_2\f$.
        ComplexType operator()(ComplexType z1, ComplexType z2) const;

        /// Add a non-resonant fermion-fermion term to this term.
        ///
        /// This operator does not check similarity of the terms!
        /// Parameters of this term are updated as follows.
        /// \li Coeff += AnotherTerm.Coeff
        /// \li Poles[i] = (Poles[i] * Weight + AnotherTerm.Poles[i] * AnotherTerm.Weight) /
        ///                (Weight + AnotherTerm.Weight)
        /// \li Weight += AnotherTerm.Weight
        /// \param[in] AnotherTerm Term to add.
        NonResonantFFTerm& operator+=(NonResonantFFTerm const& AnotherTerm);

        /// Create and commit an MPI datatype for \ref NonResonantFFTerm.
        static MPI_Datatype mpi_datatype();
    };

    /// \brief A non-resonant fermion-boson term in the Lehmann representation of \ref ThreePointSusceptibility.
    ///
    /// It is parametrized by a complex coefficient \f$C\f$, positions of real poles \f$P_1, P_{12}\f$
    /// and a coefficient \f$\xi\f$ that controls how the bosonic frequency is computed.
    /// An explicit expression for the term reads \f$\frac{C}{(z_1-P_1)(z_1 - \xi z_2 - P_{12})}\f$.
    struct NonResonantFBTerm {
        /// Coefficient \f$C\f$.
        ComplexType Coeff;

        /// Pole \f$P_1\f$.
        RealType P1 = 0;

        /// Pole \f$P_{12}\f$.
        RealType P12 = 0;

        /// Coefficient \f$\xi\f$.
        int xi = 1;

        /// Weight \f$W\f$ used in addition of terms with different poles.
        /// \see \ref operator+=()
        long Weight = 0;

        /// Hasher for non-resonant fermion-boson terms.
        struct Hash {
            /// Poles located within this energy spacing from each other produce the same hash value.
            double EnergySpacing;
            /// Constructor.
            /// \param[in] EnergySpacing Energy spacing.
            explicit Hash(double EnergySpacing = 1e-8) : EnergySpacing(EnergySpacing) {}
            /// Compute hash of a term.
            /// \param[in] t Term to compute hash for.
            std::size_t operator()(NonResonantFBTerm const& t) const {
                auto h = boost::hash<std::tuple<int, std::size_t, std::size_t>>();
                return h(std::make_tuple(t.xi,
                                         hash_binned_real(t.P1, EnergySpacing),
                                         hash_binned_real(t.P12, EnergySpacing)));
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&EnergySpacing, 1, MPI_DOUBLE, root, comm); }
        };

        /// Similarity predicate for non-resonant fermion-fermion terms.
        struct KeyEqual {
            /// Tolerance level used to compare positions of the poles.
            double Tolerance;
            /// Constructor.
            /// \param[in] Tolerance Tolerance level used to compare positions of the pole.
            explicit KeyEqual(double Tolerance = 1e-8) : Tolerance(Tolerance) {}
            /// Are terms similar?
            /// \param[in] t1 First term.
            /// \param[in] t2 Second term.
            bool operator()(NonResonantFBTerm const& t1, NonResonantFBTerm const& t2) const {
                return t2.xi == t1.xi && (std::abs(t2.P1 - t1.P1) < Tolerance) &&
                       (std::abs(t2.P12 - t1.P12) < Tolerance);
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
            bool operator()(NonResonantFBTerm const& t, std::size_t ToleranceDivisor) const {
                return std::abs(t.Coeff) < Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        NonResonantFBTerm() = default;

        /// Constructor.
        /// \param[in] Coeff Coefficient of the term \f$C\f$.
        /// \param[in] P1 Pole \f$P_1\f$.
        /// \param[in] P12 Pole \f$P_{12}\f$.
        /// \param[in] xi Coefficient \f$\xi\f$.
        inline NonResonantFBTerm(ComplexType Coeff, RealType P1, RealType P12, int xi)
            : Coeff(Coeff), P1(P1), P12(P12), xi(xi), Weight(1) {}

        /// Substitute complex frequencies \f$z_1, z_2\f$ into this term.
        /// \param[in] z1 Complex frequency \f$z_1\f$.
        /// \param[in] z2 Complex frequency \f$z_2\f$.
        ComplexType operator()(ComplexType z1, ComplexType z2) const;

        /// Add a non-resonant fermion-boson term to this term.
        ///
        /// This operator does not check similarity of the terms!
        /// Parameters of this term are updated as follows.
        /// \li Coeff += AnotherTerm.Coeff
        /// \li P = (P1 * Weight + AnotherTerm.P1 * AnotherTerm.Weight) / (Weight + AnotherTerm.Weight)
        /// \li P12 = (P12 * Weight + AnotherTerm.P12 * AnotherTerm.Weight) / (Weight + AnotherTerm.Weight)
        /// \li Weight += AnotherTerm.Weight
        /// \param[in] AnotherTerm Term to add.
        NonResonantFBTerm& operator+=(NonResonantFBTerm const& AnotherTerm);

        /// Create and commit an MPI datatype for \ref NonResonantFBTerm.
        static MPI_Datatype mpi_datatype();
    };

    /// \brief A resonant term in the Lehmann representation of \ref ThreePointSusceptibility.
    ///
    /// It is parametrized by a complex coefficient \f$C\f$, position of a real pole \f$P\f$
    /// and a coefficient \f$\xi\f$ that controls how the bosonic frequency is computed.
    /// An explicit expression for the term reads \f$\frac{C}{z_1-P}\delta_{z_1, \xi z_2}\f$.
    struct ResonantTerm {
        /// Coefficient \f$C\f$.
        ComplexType Coeff = 0;

        /// Pole \f$P\f$.
        RealType P = 0;

        /// Coefficient \f$\xi\f$.
        int xi = 1;

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
                auto h = boost::hash<std::tuple<int, std::size_t>>();
                return h(std::make_tuple(t.xi, hash_binned_real(t.P, EnergySpacing)));
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
                return t2.xi == t1.xi && (std::abs(t2.P - t1.P) < Tolerance);
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
                return std::abs(t.Coeff) < Tolerance / ToleranceDivisor;
            }
            /// Broadcast this object from a root MPI rank to all other ranks in a communicator.
            /// \param[in] comm The MPI communicator for the broadcast operation.
            /// \param[in] root Rank of the root MPI process.
            void broadcast(MPI_Comm const& comm, int root) { MPI_Bcast(&Tolerance, 1, MPI_DOUBLE, root, comm); }
        };

        ResonantTerm() = default;

        /// Constructor.
        /// \param[in] Coeff Coefficient of the term \f$C\f$.
        /// \param[in] P Pole \f$P\f$.
        /// \param[in] xi Coefficient \f$\xi\f$.
        inline ResonantTerm(ComplexType Coeff, RealType P, int xi) : Coeff(Coeff), P(P), xi(xi), Weight(1) {}

        /// Substitute complex frequencies \f$z_1, z_2\f$ into this term.
        /// \param[in] z1 Complex frequency \f$z_1\f$.
        /// \param[in] z2 Complex frequency \f$z_2\f$.
        /// \param[in] DeltaTolerance Tolerance for the resonance detection.
        ComplexType operator()(ComplexType z1, ComplexType z2, RealType DeltaTolerance = 1e-16) const;

        /// Add a resonant term to this term.
        ///
        /// This operator does not check similarity of the terms!
        /// Parameters of this term are updated as follows.
        /// \li Coeff += AnotherTerm.Coeff
        /// \li P = (P * Weight + AnotherTerm.P * AnotherTerm.Weight) / (Weight + AnotherTerm.Weight)
        /// \li Weight += AnotherTerm.Weight
        /// \param[in] AnotherTerm Term to add.
        ResonantTerm& operator+=(ResonantTerm const& AnotherTerm);

        /// Create and commit an MPI datatype for \ref ResonantTerm.
        static MPI_Datatype mpi_datatype();
    };

private:
    /// Part of the first fermionic operator.
    MonomialOperatorPart const& F1;
    /// Part of the second fermionic operator.
    MonomialOperatorPart const& F2;
    /// First multiplier of the quadratic operator \f$\hat B\f$.
    MonomialOperatorPart const& B1;
    /// Second multiplier of the quadratic operator \f$\hat B\f$.
    MonomialOperatorPart const& B2;

    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_1}\f$.
    HamiltonianPart const& Hpart1;
    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_2}\f$.
    HamiltonianPart const& Hpart2;
    /// Diagonal block of the Hamiltonian corresponding to the subspace \f${\rm S_3}\f$.
    HamiltonianPart const& Hpart3;

    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_1}\f$.
    DensityMatrixPart const& DMpart1;
    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_2}\f$.
    DensityMatrixPart const& DMpart2;
    /// Diagonal block of the many-body density matrix corresponding to the subspace \f${\rm S_3}\f$.
    DensityMatrixPart const& DMpart3;

    /// Channel
    Channel channel;

    /// Are fermionic operators in this part swapped with respect to their order in the definition?
    bool SwappedFermionOps;

    /// List of all non-resonant fermion-fermion terms contributing to this part.
    TermList<NonResonantFFTerm> NonResonantFFTerms;
    /// List of all non-resonant fermion-boson terms contributing to this part.
    TermList<NonResonantFBTerm> NonResonantFBTerms;
    /// List of all resonant terms contributing to this part.
    TermList<ResonantTerm> ResonantTerms;

    /// Adds a multi-term that has one of the following forms:
    /// \li PP channel, non-swapped fermionic operators: \f$C f(z_1, z_2)\f$;
    /// \li PP channel, swapped fermionic operators: \f$C f(z_2, z_1)\f$;
    /// \li PH/xPH channels, non-swapped fermionic operators: \f$C f(z_1, -z_2)\f$;
    /// \li PH/xPH channels, swapped fermionic operators: \f$C f(-z_2, z_1)\f$.
    ///
    /// Here, function \f$f(z_1, z_2)\f$ is defined as
    /// \f[
    /// f(z_1, z_2) = \frac{1}{z_2-(E_j-E_k)}
    ///     \left(-\frac{w_i+w_j}{z_1-(E_i-E_j)} +
    ///     \frac{(w_i-w_k)}{(z_1+z_2)-(E_i-E_k)}(1-\delta_{E_i,E _k}) +
    ///     \beta w_i \delta_{z_1+z_2,0} \delta_{E_i,E_k}\right).
    /// \f]
    ///
    /// \param[in] Coeff Common prefactor \f$C\f$.
    /// \param[in] beta Inverse temperature.
    /// \param[in] Ei The first energy level \f$E_i\f$.
    /// \param[in] Ej The second energy level \f$E_j\f$.
    /// \param[in] Ek The third energy level \f$E_k\f$.
    /// \param[in] Wi The first weight \f$w_i\f$.
    /// \param[in] Wj The second weight \f$w_j\f$.
    /// \param[in] Wk The third weight \f$w_k\f$.
    void addMultiterm(ComplexType Coeff,
                      RealType beta,
                      RealType Ei,
                      RealType Ej,
                      RealType Ek,
                      RealType Wi,
                      RealType Wj,
                      RealType Wk);

    /// A difference in energies with magnitude below this value is treated as zero.
    RealType ReduceResonanceTolerance = 1e-8;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account.
    RealType CoefficientTolerance = 1e-16;

    // compute() implementation details.
    template <bool Complex> void computeImpl();

public:
    /// Constructor.
    /// \param[in] F1 Part of the first fermionic field operator \f$\hat F_1\f$.
    /// \param[in] F2 Part of the second fermionic field operator \f$\hat F_2\f$.
    /// \param[in] B1 Part of the first multiplier of the quadratic operator \f$\hat B\f$.
    /// \param[in] B2 Part of the second multiplier of the quadratic operator \f$\hat B\f$.
    /// \param[in] Hpart1 Part of the Hamiltonian corresponding to the subspace \f${\rm S_1}\f$.
    /// \param[in] Hpart2 Part of the Hamiltonian corresponding to the subspace \f${\rm S_2}\f$.
    /// \param[in] Hpart3 Part of the Hamiltonian corresponding to the subspace \f${\rm S_3}\f$.
    /// \param[in] DMpart1 Part of the many-body density matrix corresponding to the subspace \f${\rm S_1}\f$.
    /// \param[in] DMpart2 Part of the many-body density matrix corresponding to the subspace \f${\rm S_2}\f$.
    /// \param[in] DMpart3 Part of the many-body density matrix corresponding to the subspace \f${\rm S_3}\f$.
    /// \param[in] channel Susceptibility channel.
    /// \param[in] SwappedFermionOps Are fermionic operators swapped with respect to their order in the definition?
    ThreePointSusceptibilityPart(MonomialOperatorPart const& F1,
                                 MonomialOperatorPart const& F2,
                                 MonomialOperatorPart const& B1,
                                 MonomialOperatorPart const& B2,
                                 HamiltonianPart const& Hpart1,
                                 HamiltonianPart const& Hpart2,
                                 HamiltonianPart const& Hpart3,
                                 DensityMatrixPart const& DMpart1,
                                 DensityMatrixPart const& DMpart2,
                                 DensityMatrixPart const& DMpart3,
                                 Channel channel,
                                 bool SwappedFermionOps);

    /// Compute the terms contributing to this part.
    void compute();

    /// Purge all terms.
    void clear();

    /// Substitute complex frequencies \f$z_1, z_2\f$ into this part.
    /// \param[in] z1 First frequency \f$z_1\f$.
    /// \param[in] z2 Second frequency \f$z_2\f$.
    ComplexType operator()(ComplexType z1, ComplexType z2) const;
    /// Substitute Matsubara frequencies \f$i\omega_{n_1}, i\omega_{n_2}\f$ into this part.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2) const;

    /// Return the number of resonant terms.
    std::size_t getNumResonantTerms() const { return ResonantTerms.size(); }
    /// Return the number of non-resonant fermion-fermion terms.
    std::size_t getNumNonResonantFFTerms() const { return NonResonantFFTerms.size(); }
    /// Return the number of non-resonant fermion-boson terms.
    std::size_t getNumNonResonantFBTerms() const { return NonResonantFBTerms.size(); }

    /// Are fermionic operators in this part swapped with respect to their order in the definition?
    bool getSwappedFermionOps() const { return SwappedFermionOps; }

    /// Access the list of the resonant terms.
    TermList<ResonantTerm> const& getResonantTerms() const { return ResonantTerms; }
    /// Access the list of the non-resonant fermion-fermion terms.
    TermList<NonResonantFFTerm> const& getNonResonantFFTerms() const { return NonResonantFFTerms; }
    /// Access the list of the non-resonant fermion-boson terms.
    TermList<NonResonantFBTerm> const& getNonResonantFBTerms() const { return NonResonantFBTerms; }
};

///@}

inline ComplexType ThreePointSusceptibilityPart::NonResonantFFTerm::operator()(ComplexType z1, ComplexType z2) const {
    return Coeff / ((z1 - Poles[0]) * (z2 - Poles[1]));
}

inline ComplexType ThreePointSusceptibilityPart::NonResonantFBTerm::operator()(ComplexType z1, ComplexType z2) const {
    return Coeff / ((z1 - P1) * (z1 - double(xi) * z2 - P12));
}

inline ComplexType
ThreePointSusceptibilityPart::ResonantTerm::operator()(ComplexType z1, ComplexType z2, RealType DeltaTolerance) const {
    return (std::abs(z1 - double(xi) * z2) < DeltaTolerance) ? (Coeff / (z1 - P)) : 0;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITYPART_HPP
