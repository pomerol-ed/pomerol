//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Susceptibility.hpp
/// \brief Dynamical susceptibility in the Matsubara representation.
/// \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_SUSCEPTIBILITY_HPP
#define POMEROL_INCLUDE_SUSCEPTIBILITY_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "EnsembleAverage.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "SusceptibilityPart.hpp"
#include "Thermal.hpp"

#include <complex>
#include <numeric>
#include <vector>

namespace Pomerol {

/// \defgroup Susc Dynamical susceptibilities
///@{

/// \brief Dynamical susceptibility.
///
/// This class represents a dynamical susceptibility in the Matsubara space,
/// \f[
///   \chi(i\omega_n) = \int_0^\beta Tr[\mathcal{T}_\tau \hat\rho A(\tau) B(0)] e^{i\omega_n\tau} d\tau
/// \f]
/// and its connected part,
/// \f[
///   \tilde{\chi}(i\omega_n) = \chi(i\omega_n) - \beta \langle\hat A \rangle \langle\hat B \rangle.
/// \f]
/// Here, \f$\beta\f$ is inverse temperature and \f$\langle\hat A\rangle, \langle\hat B\rangle\f$ are
/// EnsembleAverage's of boson-like monomial operators \f$\hat A, \hat B\f$.
///
/// It is actually a container class for a collection of \ref SusceptibilityPart's
/// (most of the real calculations take place in the parts).
class Susceptibility : public Thermal, public ComputableObject {

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// The Hamiltonian.
    Hamiltonian const& H;
    /// Monomial operator \f$\hat A\f$.
    MonomialOperator const& A;
    /// Monomial operator \f$\hat B\f$.
    MonomialOperator const& B;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;

    /// A flag that marks an identically vanishing susceptibility.
    bool Vanishing = true;

    /// The list of all \ref SusceptibilityPart's contributing to this susceptibility.
    std::vector<SusceptibilityPart> parts;

    /// Subtract the disconnected part \f$\langle\hat A \rangle \langle\hat B \rangle\f$?
    bool SubtractDisconnected = false;

    /// Ensemble averages \f$\langle\hat A \rangle, \langle\hat B \rangle\f$.
    ComplexType ave_A = {}, ave_B = {};

public:
    /// Constructor.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] A Monomial operator \f$\hat A\f$.
    /// \param[in] B Monomial operator \f$\hat B\f$.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    Susceptibility(StatesClassification const& S,
                   Hamiltonian const& H,
                   MonomialOperator const& A,
                   MonomialOperator const& B,
                   DensityMatrix const& DM);

    /// Copy-constructor.
    /// \param[in] Chi \ref Susceptibility object to be copied.
    Susceptibility(Susceptibility const& Chi);

    /// Select all relevant parts of \f$\hat A\f$ and \f$\hat B\f$
    /// and allocate resources for the \ref SusceptibilityPart's.
    void prepare();

    /// Actually computes the parts.
    void compute();

    /// Activate subtraction of the disconnected part \f$\langle\hat A \rangle \langle\hat B \rangle\f$.
    /// \f$\langle\hat A \rangle\f$ and \f$\langle\hat B \rangle\f$ are computed by this class.
    void subtractDisconnected();

    /// Activate subtraction of the disconnected part \f$\langle\hat A \rangle \langle\hat B \rangle\f$.
    /// Precomputed ensemble averages must be provided by the caller.
    /// \param[in] ave_A Precomputed value of \f$\langle\hat A \rangle\f$.
    /// \param[in] ave_B Precomputed value of \f$\langle\hat B \rangle\f$.
    void subtractDisconnected(ComplexType ave_A, ComplexType ave_B);

    /// Activate subtraction of the disconnected part \f$\langle\hat A \rangle \langle\hat B \rangle\f$.
    /// EnsembleAverage objects must must be provided by the caller.
    /// \param[in] EA_A Predefined EnsembleAverage class for the operator \f$\hat A\f$.
    /// \param[in] EA_B Predefined EnsembleAverage class for the operator \f$\hat B\f$.
    void subtractDisconnected(EnsembleAverage& EA_A, EnsembleAverage& EA_B);

    /// Return the susceptibility value calculated at a given Matsubara frequency.
    /// \param[in] MatsubaraNumber Index of the Matsubara frequency \f$n\f$ (\f$\omega_n=2\pi n/\beta\f$).
    ComplexType operator()(long MatsubaraNumber) const;

    /// Return the susceptibility value calculated at a given complex frequency \f$z\f$.
    /// \param[in] z The complex frequency.
    ComplexType operator()(ComplexType z) const;

    /// Return the susceptibility value calculated at a given imaginary time \f$\tau\f$.
    /// \param[in] tau Imaginary time point.
    ComplexType of_tau(RealType tau) const;

    /// Is this susceptibility identically zero?
    bool isVanishing() const { return Vanishing; }
};

///@}

inline ComplexType Susceptibility::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber));
}

inline ComplexType Susceptibility::operator()(ComplexType z) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        Value = std::accumulate(parts.begin(),
                                parts.end(),
                                ComplexType(0),
                                [z](ComplexType s, SusceptibilityPart const& p) { return s + p(z); });
    }
    if(SubtractDisconnected && std::abs(z) < 1e-15)
        Value -= ave_A * ave_B * beta; // only for n=0
    return Value;
}

inline ComplexType Susceptibility::of_tau(RealType tau) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        Value = std::accumulate(parts.begin(),
                                parts.end(),
                                ComplexType(0),
                                [tau](ComplexType s, SusceptibilityPart const& p) { return s + p.of_tau(tau); });
    }
    if(SubtractDisconnected)
        Value -= ave_A * ave_B;
    return Value;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_SUSCEPTIBILITY_HPP
