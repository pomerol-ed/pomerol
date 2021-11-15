//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/Thermal.h
** \brief Thermal object (an object which has sense only for a finite temperature).
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_THERMAL_HPP
#define POMEROL_INCLUDE_THERMAL_HPP

#include "Misc.hpp"

#include <cmath>

namespace Pomerol {

struct Thermal {
    RealType const beta;
    ComplexType const MatsubaraSpacing;

    explicit Thermal(RealType beta) : beta(beta), MatsubaraSpacing(I * M_PI / beta) {}
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_THERMAL_HPP
