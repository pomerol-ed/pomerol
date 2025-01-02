//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Vertex4.hpp
/// \brief Irreducible two-particle vertex in the Matsubara representation.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_VERTEX4_HPP
#define POMEROL_INCLUDE_VERTEX4_HPP

#include "ComputableObject.hpp"
#include "GreensFunction.hpp"
#include "MatsubaraContainers.hpp"
#include "Misc.hpp"
#include "Thermal.hpp"
#include "TwoParticleGF.hpp"

namespace Pomerol {

/// \addtogroup 2PGF
///@{

/// \brief Irreducible two-particle vertex.
///
/// Irreducible two-particle vertex part of fermions,
/// \f[
///   \Gamma_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_4}) =
///   \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_4}) -
///   \chi^0_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_4})
/// \f]
/// with the Wick part of the two-particle Green's function being
/// \f[
///   \chi^0_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_4}) =
///   \beta\delta_{\omega_{n_1}\omega_{n_4}}\delta_{\omega_{n_2}\omega_{n_3}}G_{il}(\omega_{n_1})G_{jk}(\omega_{n_2}) -
///   \beta\delta_{\omega_{n_1}\omega_{n_3}}\delta_{\omega_{n_2}\omega_{n_4}}G_{ik}(\omega_{n_1})G_{jl}(\omega_{n_2}).
/// \f]
/// \f$\beta\f$ is the inverse temperature, \f$G_{ij}\f$ is the single-particle Green's function, and
/// \f$\omega_{n_4} = \omega_{n_1}+\omega_{n_2}-\omega_{n_3}\f$.
class Vertex4 : public Thermal, public ComputableObject {

    /// Fermionic two-particle Matsubara Green's function \f$\chi_{ijkl}\f$.
    TwoParticleGF const& Chi;
    /// Fermionic single-particle Matsubara Green's function \f$G_{ik}\f$.
    GreensFunction const& G13;
    /// Fermionic single-particle Matsubara Green's function \f$G_{jl}\f$.
    GreensFunction const& G24;
    /// Fermionic single-particle Matsubara Green's function \f$G_{il}\f$.
    GreensFunction const& G14;
    /// Fermionic single-particle Matsubara Green's function \f$G_{jk}\f$.
    GreensFunction const& G23;

    /// Storage for precomputed values.
    mutable MatsubaraContainer4<Vertex4> Storage;

    friend class MatsubaraContainer4<Vertex4>;

public:
    /// Constructor.
    /// \param[in] Chi Fermionic two-particle Matsubara Green's function \f$\chi_{ijkl}\f$
    /// \param[in] G13 Fermionic single-particle Matsubara Green's function \f$G_{ik}\f$.
    /// \param[in] G24 Fermionic single-particle Matsubara Green's function \f$G_{jl}\f$.
    /// \param[in] G14 Fermionic single-particle Matsubara Green's function \f$G_{il}\f$.
    /// \param[in] G23 Fermionic single-particle Matsubara Green's function \f$G_{jk}\f$.
    Vertex4(TwoParticleGF const& Chi,
            GreensFunction const& G13,
            GreensFunction const& G24,
            GreensFunction const& G14,
            GreensFunction const& G23);

    /// Populate the internal cache of precomputed values.
    /// \param[in] NumberOfMatsubaras Number of positive fermionic Matsubara frequencies
    ///            \f$\omega_{n_1}\f$ and \f$\omega_{n_2}\f$ for which values are precomputed and stored.
    void compute(long NumberOfMatsubaras = 0);

    /// Return the value of the vertex calculated a given Matsubara frequency triplet.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    /// \param[in] MatsubaraNumber3 Index of the third Matsubara frequency
    ///                             \f$n_3\f$ (\f$\omega_{n_3}=\pi(2n_3+1)/\beta\f$).
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /// Return the value of vertex calculated a given Matsubara frequency triplet.
    /// This method ignores the internal cache of precomputed values.
    /// \param[in] MatsubaraNumber1 Index of the first Matsubara frequency
    ///                             \f$n_1\f$ (\f$\omega_{n_1}=\pi(2n_1+1)/\beta\f$).
    /// \param[in] MatsubaraNumber2 Index of the second Matsubara frequency
    ///                             \f$n_2\f$ (\f$\omega_{n_2}=\pi(2n_2+1)/\beta\f$).
    /// \param[in] MatsubaraNumber3 Index of the third Matsubara frequency
    ///                             \f$n_3\f$ (\f$\omega_{n_3}=\pi(2n_3+1)/\beta\f$).
    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_VERTEX4_HPP
