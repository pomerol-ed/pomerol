/** \file Index.h
**  \brief Types related to indices and their combinations.
**
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef POMEROL_INCLUDE_POMEROL_INDEX_H
#define POMEROL_INCLUDE_POMEROL_INDEX_H

#include"Misc.h"

#include <ostream>
#include <tuple>

namespace Pomerol {

/** A structure to handle a combination of 2 Particle Indices. */
struct IndexCombination2
{
    /** Actual indices */
    ParticleIndex Index1, Index2;

    /** Constructor
     * \param[in] Index1 - Index of a 1st operator (C)
     * \param[in] Index2 - Index of a 2nd operator (C^+)
    */
    IndexCombination2(ParticleIndex Index1, ParticleIndex Index2) :
        Index1(Index1), Index2(Index2) {}

    /** Operator < - comparison method for IndexCombination */
    bool operator< (const IndexCombination2& rhs) const
    {
        return std::make_tuple(Index1, Index2) <
               std::make_tuple(rhs.Index1, rhs.Index2);
    }
    /** Operator == */
    bool operator==(const IndexCombination2& rhs) const
    {
        return Index1 == rhs.Index1 && Index2 == rhs.Index2;
    }
    /** Operator != */
    bool operator!=(const IndexCombination2& rhs) const
    {
        return !operator==(rhs);
    }

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& os, const IndexCombination2& ic)
    {
        return os << "(" << ic.Index1 << ic.Index2 << ")";
    }
};

/** A structure to handle a combination of 4 Particle Indices. */
struct IndexCombination4
{
    /** Actual indices */
    ParticleIndex Index1, Index2, Index3, Index4;

    /** Constructor
     * \param[in] Index1 - Index of a 1st operator (C)
     * \param[in] Index2 - Index of a 2nd operator (C)
     * \param[in] Index3 - Index of a 3nd operator (C^+)
     * \param[in] Index4 - Index of a 4nd operator (C^+)
    */
    IndexCombination4(ParticleIndex Index1, ParticleIndex Index2,
                     ParticleIndex Index3, ParticleIndex Index4) :
        Index1(Index1), Index2(Index2), Index3(Index3), Index4(Index4) {}

    /** Operator < - comparison method for IndexCombination */
    bool operator< (const IndexCombination4& rhs) const
    {
        return std::make_tuple(Index1, Index2, Index3, Index4) <
               std::make_tuple(rhs.Index1, rhs.Index2, rhs.Index3, rhs.Index4);
    }
    /** Operator == */
    bool operator==(const IndexCombination4& rhs) const
    {
        return Index1 == rhs.Index1 &&
               Index2 == rhs.Index2 &&
               Index3 == rhs.Index3 &&
               Index4 == rhs.Index4;
    }
    /** Operator != */
    bool operator!=(const IndexCombination4& rhs) const
    {
        return !operator==(rhs);
    }

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& os, const IndexCombination4& ic)
    {
        return os << "(" << ic.Index1 << ic.Index2 << ic.Index3 << ic.Index4 << ")";
    }
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_INDEX_H
