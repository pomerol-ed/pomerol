//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Thermal.hpp
/// \brief Objects, whose definition depends on the temperature.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_THERMAL_HPP
#define POMEROL_INCLUDE_THERMAL_HPP

#include "Misc.hpp"

#include <cmath>

namespace Pomerol {

/// \addtogroup Misc
///@{

/// Base class for objects whose definition depends on the temperature.
struct Thermal {
    /// Inverse temperature \f$\beta\f$.
    RealType const beta;
    /// Spacing between (imaginary) Matsubara frequencies, \f$i\pi/\beta\f$.
    ComplexType const MatsubaraSpacing;

    /// Construct a thermal object for a given inverse temperature.
    /// \param[in] beta Inverse temperature \f$\beta\f$
    explicit Thermal(RealType beta) : beta(beta), MatsubaraSpacing(I * M_PI / beta) {}
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THERMAL_HPP
