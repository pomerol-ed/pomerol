//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.
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

// Catch2 matcher class that checks proximity of two complex numbers
class IsCloseToMatcher : public Catch::MatcherBase<std::complex<double>> {
    std::complex<double> ref;
    double tol;

public:

    IsCloseToMatcher(std::complex<double> ref, double tol) :
        ref(std::move(ref)), tol(tol)
    {}

    bool match(std::complex<double> const& x) const override {
        return std::abs(x - ref) <= tol;
    }

    std::string describe() const override {
        std::ostringstream ss;
        ss << "is close to " << ref << " (tol = " << tol << ")";
        return ss.str();
    }
};

inline IsCloseToMatcher
IsCloseTo(std::complex<double> const& ref, double tol = 1e-10) {
    return IsCloseToMatcher(ref, tol);
}

#endif // #ifndef POMEROL_TEST_CATCH2_CATCH_POMEROL_HPP
