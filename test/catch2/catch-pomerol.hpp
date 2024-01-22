//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/catch2/catch-pomerol.hpp
/// \brief Supporting types and functions for Catch2
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_TEST_CATCH2_CATCH_POMEROL_HPP
#define POMEROL_TEST_CATCH2_CATCH_POMEROL_HPP

// Undefine INFO() from Misc.hpp as Catch2 defines its own INFO() macro,
// which is more appropriate for unit testing.
#ifdef INFO
#undef INFO
#endif

#include "catch.hpp"

#include <cmath>
#include <complex>
#include <sstream>
#include <utility>

/// Catch2 matcher class that checks proximity of two complex numbers
class IsCloseToMatcher : public Catch::MatcherBase<std::complex<double>> {
    std::complex<double> ref;
    double tol;

public:
    IsCloseToMatcher(std::complex<double> ref, double tol) : ref(std::move(ref)), tol(tol) {}

    bool match(std::complex<double> const& x) const override { return std::abs(x - ref) <= tol; }

    std::string describe() const override {
        std::ostringstream ss;
        ss << "is close to " << ref << " (tol = " << tol << ")";
        return ss.str();
    }
};

/// Factory function for \ref IsCloseToMatcher matchers.
/// \param[in] ref Reference complex value to compare to.
/// \param[in] tol Maximum absolute deviation from the reference value.
inline IsCloseToMatcher IsCloseTo(std::complex<double> const& ref, double tol = 1e-10) {
    return IsCloseToMatcher(ref, tol);
}

#endif // #ifndef POMEROL_TEST_CATCH2_CATCH_POMEROL_HPP
