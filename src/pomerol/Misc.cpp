//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/Misc.cpp
/// \brief Implementation of functions and methods from Misc.hpp
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/Misc.hpp"

#include <cmath>
#include <cstdint>
#include <functional>

namespace Pomerol {

//////////////////
// Permutation3 //
//////////////////
bool Permutation3::operator==(Permutation3 const& rhs) const {
    return (sign == rhs.sign && perm[0] == rhs.perm[0] && perm[1] == rhs.perm[1]);
}

bool Permutation3::operator!=(Permutation3 const& rhs) const {
    return !(*this == rhs);
}

std::ostream& operator<<(std::ostream& out, Permutation3 const& p) {
    if(p.sign == -1)
        out << "-";
    return out << p.perm[0] + 1 << p.perm[1] + 1 << p.perm[2] + 1;
}

std::array<Permutation3, 6> const permutations3 = {
    {{{0, 1, 2}, 1}, {{0, 2, 1}, -1}, {{1, 0, 2}, -1}, {{1, 2, 0}, 1}, {{2, 0, 1}, 1}, {{2, 1, 0}, -1}}};

//////////////////
// Permutation4 //
//////////////////
bool Permutation4::operator==(Permutation4 const& rhs) const {
    return (sign == rhs.sign && perm[0] == rhs.perm[0] && perm[1] == rhs.perm[1] && perm[2] == rhs.perm[2]);
}

bool Permutation4::operator!=(Permutation4 const& rhs) const {
    return !(*this == rhs);
}

std::ostream& operator<<(std::ostream& out, Permutation4 const& p) {
    if(p.sign == -1)
        out << "-";
    return out << p.perm[0] + 1 << p.perm[1] + 1 << p.perm[2] + 1 << p.perm[3] + 1;
}

std::array<Permutation4, 24> const permutations4 = {
    {{{0, 1, 2, 3}, 1},  {{0, 1, 3, 2}, -1}, {{0, 2, 1, 3}, -1}, {{0, 2, 3, 1}, 1},  {{0, 3, 1, 2}, 1},
     {{0, 3, 2, 1}, -1}, {{1, 0, 2, 3}, -1}, {{1, 0, 3, 2}, 1},  {{1, 2, 0, 3}, 1},  {{1, 2, 3, 0}, -1},
     {{1, 3, 0, 2}, -1}, {{1, 3, 2, 0}, 1},  {{2, 0, 1, 3}, 1},  {{2, 0, 3, 1}, -1}, {{2, 1, 0, 3}, -1},
     {{2, 1, 3, 0}, 1},  {{2, 3, 0, 1}, 1},  {{2, 3, 1, 0}, -1}, {{3, 0, 1, 2}, -1}, {{3, 0, 2, 1}, 1},
     {{3, 1, 0, 2}, 1},  {{3, 1, 2, 0}, -1}, {{3, 2, 0, 1}, -1}, {{3, 2, 1, 0}, 1}}};

/////////////
// Channel //
/////////////
std::ostream& operator<<(std::ostream& os, Channel channel) {
    switch(channel) {
    case Channel::PP: return os << "PP";
    case Channel::PH: return os << "PH";
    case Channel::xPH: return os << "xPH";
    default: return os;
    }
}

////////////////////////
// hash_binned_real() //
////////////////////////
std::size_t hash_binned_real(double x, double bin_size) {
    auto h = std::hash<std::int64_t>();
    return h(std::round(x / bin_size));
}

} // namespace Pomerol
