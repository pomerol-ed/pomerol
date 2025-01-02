//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/ThreePointSusceptibility.hpp
/// \brief 3-point susceptibility in the Matsubara representation.
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITY_HPP
#define POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITY_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"
#include "ThreePointSusceptibilityPart.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <cstddef>
#include <numeric>
#include <tuple>
#include <vector>

namespace Pomerol {

/// \defgroup 3PSusc Three-point fermion-boson susceptibility
///@{

/// Duplet of complex frequencies.
using FreqTuple2 = std::tuple<ComplexType, ComplexType>;
/// List of complex frequency duplets.
using FreqVec2 = std::vector<FreqTuple2>;

/// \brief 3-point fermion-boson susceptibility in the Matsubara representation.
///
/// The susceptibility can be defined in one of the following three channels.
/// \li Particle-particle channel: \f[\chi^{(3)}_{pp}(\omega_{n_1},\omega_{n_2}) =
/// \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{-i\omega_{n_2}\tau_2}
/// Tr[\mathcal{T}_\tau \hat\rho c^\dagger_1(\tau_1) c_2(0^+) c^\dagger_3(\tau_2) c_4(0)]\f]
///
/// \li Particle-hole channel: \f[\chi^{(3)}_{ph}(\omega_{n_1},\omega_{n_2}) =
/// \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{i\omega_{n_2}\tau_2}
/// Tr[\mathcal{T}_\tau \hat\rho c^\dagger_1(\tau_1) c_2(\tau_2) c^\dagger_3(0^+) c_4(0)]\f]
///
/// \li Crossed particle-hole channel: \f[\chi^{(3)}_{\bar{ph}}(\omega_{n_1},\omega_{n_2}) =
/// \int_0^\beta d\tau_1 d\tau_2 e^{-i\omega_{n_1}\tau_1} e^{i\omega_{n_2}\tau_2}
/// Tr[\mathcal{T}_\tau \hat\rho c^\dagger_1(\tau_1) c_2(0) c^\dagger_3(0^+) c_4(\tau_2)]\f]
///
/// These susceptibilities can be interpreted as 3-point correlators of two fermionic
/// operators \f$\hat F_1(\omega_{n_1}), \hat F_2(\omega_{n_2})\f$ and one quadratic operator \f$\hat B\f$.
/// \li PP channel: \f$\hat F_1 = c^\dagger_1, \hat F_2 = c^\dagger_3, \hat B = \Delta_{24} = c_2 c_4 \f$;
/// \li PH channel: \f$\hat F_1 = c^\dagger_1, \hat F_2 = c_2, \hat B = n_{34} = c^\dagger_3 c_4 \f$;
/// \li xPH channel: \f$\hat F_1 = c^\dagger_1, \hat F_2 = c_4, \hat B = -n_{32} = -c^\dagger_3 c_2 \f$.
///
/// It is actually a container class for a collection of \ref ThreePointSusceptibilityPart's
/// (most of the real calculations take place in the parts).
class ThreePointSusceptibility : public Thermal, public ComputableObject {

    friend class ThreePointSusceptibilityContainer;

    /// Channel
    Channel channel;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// The Hamiltonian.
    Hamiltonian const& H;
    /// The creation operator \f$c^\dagger_1\f$.
    CreationOperator const& CX1;
    /// The annihilation operator \f$c_2\f$.
    AnnihilationOperator const& C2;
    /// The creation operator \f$c^\dagger_3\f$.
    CreationOperator const& CX3;
    /// The annihilation operator \f$c_4\f$.
    AnnihilationOperator const& C4;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;

    /// The list of all \ref ThreePointSusceptibilityPart's contributing to this susceptibility.
    std::vector<ThreePointSusceptibilityPart> parts;

protected:
    /// A flag that marks an identically vanishing susceptibility.
    bool Vanishing = true;

public:
    /// A difference in energies with magnitude below this value is treated as zero.
    RealType ReduceResonanceTolerance = 1e-8;
    /// Minimal magnitude of the coefficient of a term for it to be taken into account.
    RealType CoefficientTolerance = 1e-16;

    /// Constructor
    /// \param[in] channel Channel of the 3-point susceptibility.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] CX1 The creation operator \f$c^\dagger_1\f$.
    /// \param[in] C2 The annihilation operator \f$c_2\f$.
    /// \param[in] CX3 The creation operator \f$c^\dagger_3\f$.
    /// \param[in] C4 The annihilation operator \f$c_4\f$.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    ThreePointSusceptibility(Channel channel,
                             StatesClassification const& S,
                             Hamiltonian const& H,
                             CreationOperator const& CX1,
                             AnnihilationOperator const& C2,
                             CreationOperator const& CX3,
                             AnnihilationOperator const& C4,
                             DensityMatrix const& DM);

    /// Select relevant parts of \f$c^\dagger_1, c_2, c^\dagger_3, c_4\f$ and allocate resources for the parts.
    void prepare();

    /// Compute the parts in parallel and fill the internal cache of precomputed values.
    /// \param[in] clear If true, computed \ref ThreePointSusceptibilityPart's will be destroyed immediately after
    ///                  filling the precomputed value cache.
    /// \param[in] freqs List of frequency duplets \f$(\omega_{n_1},\omega_{n_2})\f$.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \return A list of precomputed values.
    /// \pre \ref prepare() has been called.
    std::vector<ComplexType>
    compute(bool clear = false, FreqVec2 const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);

    /// Returns the single particle index of one of the operators \f$c^\dagger_1, c_2, c^\dagger_3, c_4\f$.
    /// \param[in] Position Position of the requested operator, 0--3.
    ParticleIndex getIndex(std::size_t Position) const;

    /// Return the value of the 3-point susceptibility calculated at a given complex frequency duplet.
    /// This method ignores the precomputed value cache.
    /// \param[in] z1 First frequency \f$z_1\f$.
    /// \param[in] z2 Second frequency \f$z_2\f$.
    ComplexType operator()(ComplexType z1, ComplexType z2) const;
    /// Return the value of the 3-point susceptibility calculated a given Matsubara frequency duplet.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2) const;

    /// Is this susceptibility identically zero?
    bool isVanishing() const { return Vanishing; }

private:
    /// Select operators F1, F2, B1 and B2 depending on the selected channel
    CreationOperator const& getF1() const;
    MonomialOperator const& getF2() const;
    MonomialOperator const& getB1() const;
    MonomialOperator const& getB2() const;
};

///@}

inline ComplexType ThreePointSusceptibility::operator()(ComplexType z1, ComplexType z2) const {
    if(Vanishing)
        return 0;
    else {
        return std::accumulate(
            parts.begin(),
            parts.end(),
            ComplexType(0),
            [z1, z2](ComplexType s, ThreePointSusceptibilityPart const& p) { return s + p(z1, z2); });
    }
}

inline ComplexType ThreePointSusceptibility::operator()(long MatsubaraNumber1, long MatsubaraNumber2) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber1 + 1),
                   MatsubaraSpacing * RealType(2 * MatsubaraNumber2 + 1));
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THREEPOINTSUSCEPTIBILITY_HPP
