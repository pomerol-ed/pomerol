//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/TwoParticleGF.hpp
/// \brief Fermionic two-particle Matsubara Green's function.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_TWOPARTICLEGF_HPP
#define POMEROL_INCLUDE_TWOPARTICLEGF_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"
#include "TwoParticleGFPart.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <cstddef>
#include <numeric>
#include <tuple>
#include <vector>

namespace Pomerol {

/// \defgroup 2PGF Two-particle Green's functions of fermions
///@{

/// Triplet of complex frequencies.
using FreqTuple3 = std::tuple<ComplexType, ComplexType, ComplexType>;
/// List of complex frequency triplets.
using FreqVec3 = std::vector<FreqTuple3>;

#ifndef DOXYGEN_SKIP
using FreqTuple POMEROL_DEPRECATED = FreqTuple3;
using FreqVec POMEROL_DEPRECATED = FreqVec3;
#endif

/// \brief Fermionic two-particle Matsubara Green's function.
///
/// \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
///   \int_0^\beta
///   Tr[\mathcal{T}_\tau \hat\rho c_i(\tau_1)c_j(\tau_2)c^\dagger_k(\tau_3)c^\dagger_l(0)]
///   e^{i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3}
///   d\tau_1 d\tau_2 d\tau_3.
/// \f]
/// It is actually a container class for a collection of \ref TwoParticleGFPart's
/// (most of the real calculations take place in the parts).
class TwoParticleGF : public Thermal, public ComputableObject {

    friend class TwoParticleGFContainer;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// The Hamiltonian.
    Hamiltonian const& H;
    /// The annihilation operator \f$c_i\f$.
    AnnihilationOperator const& C1;
    /// The annihilation operator \f$c_j\f$.
    AnnihilationOperator const& C2;
    /// The creation operator \f$c^\dagger_k\f$
    CreationOperator const& CX3;
    /// The creation operator \f$c^\dagger_l\f$
    CreationOperator const& CX4;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;

    /// The list of all \ref TwoParticleGFPart's contributing to this GF.
    std::vector<TwoParticleGFPart> parts;

protected:
    /// A flag that marks an identically vanishing Green's function.
    bool Vanishing = true;

    /// Extract the operator part standing at a specified position in a given permutation of the list
    /// \f$\{c_i,c_j,c^\dagger_k,c^\dagger_l\}\f$.
    /// \param[in] PermutationNumber Serial number of the permutation within \ref permutations3.
    /// \param[in] OperatorPosition Position of the operator, 0--3.
    /// \param[in] LeftIndex The left invariant subspace index referring to the requested part.
    MonomialOperatorPart const&
    OperatorPartAtPosition(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber LeftIndex) const;

    /// Choose the operator standing at a specified position in a given permutation of the list
    /// \f$\{c_i,c_j,c^\dagger_k,c^\dagger_l\}\f$ and
    /// return its left invariant subspace index corresponding to a given right subspace index.
    /// Return \ref INVALID_BLOCK_NUMBER if the operator does not have such a (non-zero) block.
    /// \param[in] PermutationNumber Serial number of the permutation within \ref permutations3.
    /// \param[in] OperatorPosition Position of the operator, 0--3.
    /// \param[in] RightIndex The right invariant subspace index.
    BlockNumber getLeftIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber RightIndex) const;
    /// Choose the operator standing at a specified position in a given permutation of the list
    /// \f$\{c_i,c_j,c^\dagger_k,c^\dagger_l\}\f$ and
    /// return its right invariant subspace index corresponding to a given left subspace index.
    /// Return \ref INVALID_BLOCK_NUMBER if the operator does not have such a (non-zero) block.
    /// \param[in] PermutationNumber Serial number of the permutation within \ref permutations3.
    /// \param[in] OperatorPosition Position of the operator, 0--3.
    /// \param[in] LeftIndex The left invariant subspace index.
    BlockNumber getRightIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber LeftIndex) const;

public:
    /// Lehmann representation: Maximal distance between energy poles to be consider coinciding.
    RealType PoleResolution = 1e-8;
    /// Lehmann representation: Maximal magnitude of a term coefficient to be considered negligible.
    RealType CoefficientTolerance = 1e-16;

    /// Constructor.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] C1 The annihilation operator \f$c_i\f$.
    /// \param[in] C2 The annihilation operator \f$c_j\f$.
    /// \param[in] CX3 The creation operator \f$c^\dagger_k\f$.
    /// \param[in] CX4 The creation operator \f$c^\dagger_l\f$.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    TwoParticleGF(StatesClassification const& S,
                  Hamiltonian const& H,
                  AnnihilationOperator const& C1,
                  AnnihilationOperator const& C2,
                  CreationOperator const& CX3,
                  CreationOperator const& CX4,
                  DensityMatrix const& DM);

    /// Choose relevant parts of \f$c_i,c_j,c^\dagger_k,c^\dagger_l\f$ and allocate resources for the parts.
    void prepare();

    /// Compute the parts in parallel and fill the internal cache of precomputed values.
    /// \param[in] clear If true, computed \ref TwoParticleGFPart's will be destroyed immediately after
    ///                  filling the precomputed value cache.
    /// \param[in] freqs List of frequency triplets \f$(\omega_{n_1},\omega_{n_2},\omega_{n_3})\f$.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \return A list of precomputed values.
    /// \pre \ref prepare() has been called.
    std::vector<ComplexType>
    compute(bool clear = false, FreqVec3 const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);

    /// Returns the single particle index of one of the operators \f$c_i,c_j,c^\dagger_k,c^\dagger_l\f$.
    /// \param[in] Position Position of the requested operator, 0--3.
    ParticleIndex getIndex(std::size_t Position) const;

    /// Return the value of the two-particle Green's function calculated at a given complex frequency triplet.
    /// This method ignores the precomputed value cache.
    /// \param[in] z1 First frequency \f$z_1\f$.
    /// \param[in] z2 Second frequency \f$z_2\f$.
    /// \param[in] z3 Third frequency \f$z_3\f$.
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;
    /// Return the value of the two-particle Green's function calculated a given Matsubara frequency triplet.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    /// \param[in] MatsubaraNumber3 Index of the third Matsubara frequency
    ///                             \f$n_3\f$ (\f$\omega_{n_3}=\pi(2n_3+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /// Is this Green's function identically zero?
    bool isVanishing() const { return Vanishing; }
};

///@}

inline ComplexType TwoParticleGF::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const {
    if(Vanishing)
        return 0;
    else {
        return std::accumulate(parts.begin(),
                               parts.end(),
                               ComplexType(0),
                               [z1, z2, z3](ComplexType s, TwoParticleGFPart const& p) { return s + p(z1, z2, z3); });
    }
}

inline ComplexType
TwoParticleGF::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber1 + 1),
                   MatsubaraSpacing * RealType(2 * MatsubaraNumber2 + 1),
                   MatsubaraSpacing * RealType(2 * MatsubaraNumber3 + 1));
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGF_HPP
