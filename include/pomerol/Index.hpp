//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Index.hpp
/// \brief Combinations of single-particle indices.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_INDEX_HPP
#define POMEROL_INCLUDE_POMEROL_INDEX_HPP

#include "Misc.hpp"

#include <ostream>
#include <tuple>

namespace Pomerol {

/// \defgroup Misc Miscellaneous
///@{

/// A tuple-like combination of two single-particle indices.
struct IndexCombination2 {
    /// First single-particle index.
    ParticleIndex Index1;
    /// Second single-particle index.
    ParticleIndex Index2;

    /// Construct from two single-particle indices.
    /// \param[in] Index1 First single-particle index.
    /// \param[in] Index2 Second single-particle index.
    IndexCombination2(ParticleIndex Index1, ParticleIndex Index2) : Index1(Index1), Index2(Index2) {}

    /// Lexicographical less comparison operator for pairs of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator<(IndexCombination2 const& rhs) const {
        return std::make_tuple(Index1, Index2) < std::make_tuple(rhs.Index1, rhs.Index2);
    }
    /// Equal comparison operator for pairs of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator==(IndexCombination2 const& rhs) const { return Index1 == rhs.Index1 && Index2 == rhs.Index2; }
    /// Unequal comparison operator for pairs of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator!=(IndexCombination2 const& rhs) const { return !operator==(rhs); }

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] ic Index combination to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, IndexCombination2 const& ic) {
        return os << "(" << ic.Index1 << ic.Index2 << ")";
    }
};

/// A tuple-like combination of four single-particle indices.
struct IndexCombination4 {
    /// First single-particle index.
    ParticleIndex Index1;
    /// Second single-particle index.
    ParticleIndex Index2;
    /// Third single-particle index.
    ParticleIndex Index3;
    /// Fourth single-particle index.
    ParticleIndex Index4;

    /// Construct from four single-particle indices.
    /// \param[in] Index1 First single-particle index.
    /// \param[in] Index2 Second single-particle index.
    /// \param[in] Index3 Third single-particle index.
    /// \param[in] Index4 Fourth single-particle index.
    IndexCombination4(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4)
        : Index1(Index1), Index2(Index2), Index3(Index3), Index4(Index4) {}

    /// Lexicographical less comparison operator for 4-tuples of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator<(IndexCombination4 const& rhs) const {
        return std::make_tuple(Index1, Index2, Index3, Index4) <
               std::make_tuple(rhs.Index1, rhs.Index2, rhs.Index3, rhs.Index4);
    }
    /// Equal comparison operator for 4-tuples of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator==(IndexCombination4 const& rhs) const {
        return Index1 == rhs.Index1 && Index2 == rhs.Index2 && Index3 == rhs.Index3 && Index4 == rhs.Index4;
    }
    /// Unequal comparison operator for 4-tuples of indices.
    /// \param[in] rhs Index combination object to compare to.
    bool operator!=(IndexCombination4 const& rhs) const { return !operator==(rhs); }

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] ic Index combination to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, IndexCombination4 const& ic) {
        return os << "(" << ic.Index1 << ic.Index2 << ic.Index3 << ic.Index4 << ")";
    }
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEX_HPP
