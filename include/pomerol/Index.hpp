//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file Index.h
**  \brief Types related to indices and their combinations.
**
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef POMEROL_INCLUDE_POMEROL_INDEX_H
#define POMEROL_INCLUDE_POMEROL_INDEX_H

#include "Misc.hpp"

#include <ostream>
#include <tuple>

namespace Pomerol {

/** A structure to handle a combination of 2 Particle Indices. */
struct IndexCombination2 {
    /** Actual indices */
    ParticleIndex Index1, Index2;

    /** Constructor
     * \param[in] Index1 - Index of a 1st operator (C)
     * \param[in] Index2 - Index of a 2nd operator (C^+)
    */
    IndexCombination2(ParticleIndex Index1, ParticleIndex Index2) : Index1(Index1), Index2(Index2) {}

    /** Operator < - comparison method for IndexCombination */
    bool operator<(IndexCombination2 const& rhs) const {
        return std::make_tuple(Index1, Index2) < std::make_tuple(rhs.Index1, rhs.Index2);
    }
    /** Operator == */
    bool operator==(IndexCombination2 const& rhs) const { return Index1 == rhs.Index1 && Index2 == rhs.Index2; }
    /** Operator != */
    bool operator!=(IndexCombination2 const& rhs) const { return !operator==(rhs); }

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& os, IndexCombination2 const& ic) {
        return os << "(" << ic.Index1 << ic.Index2 << ")";
    }
};

/** A structure to handle a combination of 4 Particle Indices. */
struct IndexCombination4 {
    /** Actual indices */
    ParticleIndex Index1, Index2, Index3, Index4;

    /** Constructor
     * \param[in] Index1 - Index of a 1st operator (C)
     * \param[in] Index2 - Index of a 2nd operator (C)
     * \param[in] Index3 - Index of a 3nd operator (C^+)
     * \param[in] Index4 - Index of a 4nd operator (C^+)
    */
    IndexCombination4(ParticleIndex Index1, ParticleIndex Index2, ParticleIndex Index3, ParticleIndex Index4)
        : Index1(Index1), Index2(Index2), Index3(Index3), Index4(Index4) {}

    /** Operator < - comparison method for IndexCombination */
    bool operator<(IndexCombination4 const& rhs) const {
        return std::make_tuple(Index1, Index2, Index3, Index4) <
               std::make_tuple(rhs.Index1, rhs.Index2, rhs.Index3, rhs.Index4);
    }
    /** Operator == */
    bool operator==(IndexCombination4 const& rhs) const {
        return Index1 == rhs.Index1 && Index2 == rhs.Index2 && Index3 == rhs.Index3 && Index4 == rhs.Index4;
    }
    /** Operator != */
    bool operator!=(IndexCombination4 const& rhs) const { return !operator==(rhs); }

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& os, IndexCombination4 const& ic) {
        return os << "(" << ic.Index1 << ic.Index2 << ic.Index3 << ic.Index4 << ")";
    }
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEX_H
