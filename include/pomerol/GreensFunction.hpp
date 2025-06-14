//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/GreensFunction.hpp
/// \brief Fermionic single-particle Matsubara Green's function.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_GREENSFUNCTION_HPP
#define POMEROL_INCLUDE_POMEROL_GREENSFUNCTION_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "GreensFunctionPart.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <cstddef>
#include <vector>

namespace Pomerol {

/// \defgroup GF Single-particle Green's functions of fermions
///@{

/// \brief Fermionic single-particle Matsubara Green's function.
///
/// This class gives access to the GF values both in imaginary time representation,
/// \f[
///  G(\tau) = -Tr[\mathcal{T}_\tau \hat\rho c(\tau) c^\dagger(0)]
/// \f]
/// and at the imaginary Matsubara frequencies \f$\omega_n = \pi(2n+1)/\beta\f$,
/// \f[
///  G(i\omega_n) = \int_0^\beta d\tau e^{i\omega_n\tau} G(\tau).
/// \f]
/// It is actually a container class for a collection of \ref GreensFunctionPart's
/// (most of the real calculations take place in the parts).
class GreensFunction : public Thermal, public ComputableObject {

    friend class GFContainer;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// The Hamiltonian.
    Hamiltonian const& H;
    /// The annihilation operator \f$c\f$.
    AnnihilationOperator const& C;
    /// The creation operator \f$c^\dagger\f$.
    CreationOperator const& CX;
    /// Many-body density matrix \f$\hat\rho\f$.
    DensityMatrix const& DM;

    /// A flag that marks an identically vanishing Green's function.
    bool Vanishing = true;

    /// The list of all \ref GreensFunctionPart's contributing to this GF.
    std::vector<GreensFunctionPart> parts;

public:
    /// Matrix elements with magnitudes equal to or below this value are treated as negligible.
    RealType MatrixElementTolerance = 1e-8;

    /// Constructor.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian.
    /// \param[in] C The annihilation operator \f$c\f$.
    /// \param[in] CX The creation operator \f$c^\dagger\f$.
    /// \param[in] DM Many-body density matrix \f$\hat\rho\f$.
    GreensFunction(StatesClassification const& S,
                   Hamiltonian const& H,
                   AnnihilationOperator const& C,
                   CreationOperator const& CX,
                   DensityMatrix const& DM);

    /// Copy-constructor.
    /// \param[in] GF \ref GreensFunction object to be copied.
    GreensFunction(GreensFunction const& GF);

    /// Select all relevant parts of \f$c\f$ and \f$c^\dagger\f$
    /// and allocate resources for the \ref GreensFunctionPart's.
    void prepare();

    /// Actually computes the parts.
    void compute();

    /// Returns the single-particle index of either \f$c\f$ or \f$c^\dagger\f$.
    /// \param[in] Position Select \f$c\f$ for Position==0 and \f$c^\dagger\f$ for Position==1.
    ///            Throws for other values of this argument.
    ParticleIndex getIndex(std::size_t Position) const;

    /// Return the GF value calculated at a given Matsubara frequency.
    /// \param[in] MatsubaraNumber Index of the Matsubara frequency \f$n\f$ (\f$\omega_n=\pi(2n+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber) const;

    /// Return the GF value calculated at a given complex frequency \f$z\f$.
    /// \param[in] z The complex frequency.
    ComplexType operator()(ComplexType z) const;

    /// Return the GF value calculated at a given imaginary time \f$\tau\f$.
    /// \param[in] tau Imaginary time point.
    ComplexType of_tau(RealType tau) const;

    /// Is this Green's function identically zero?
    bool isVanishing() const { return Vanishing; }
};

///@}

inline ComplexType GreensFunction::operator()(long MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber + 1));
}

inline ComplexType GreensFunction::operator()(ComplexType z) const {
    if(Vanishing)
        return 0;
    else {
        return std::accumulate(parts.begin(),
                               parts.end(),
                               ComplexType(0),
                               [z](ComplexType s, GreensFunctionPart const& p) { return s + p(z); });
    }
}

inline ComplexType GreensFunction::of_tau(RealType tau) const {
    if(Vanishing)
        return 0;
    else {
        return std::accumulate(parts.begin(),
                               parts.end(),
                               ComplexType(0),
                               [tau](ComplexType s, GreensFunctionPart const& p) { return s + p.of_tau(tau); });
    }
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_GREENSFUNCTION_HPP
