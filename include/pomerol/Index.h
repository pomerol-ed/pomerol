/** \file Index.h
**  \brief Types related to indices and their combinations.
**
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef __INCLUDE_INDEX_H
#define __INCLUDE_INDEX_H

#include"Misc.h"

#include <ostream>

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
    IndexCombination2(ParticleIndex Index1, ParticleIndex Index2);

    /** Operator < - comparison method for IndexCombination */
    bool operator< (const IndexCombination2& rhs) const ;
    /** Operator == */
    bool operator==(const IndexCombination2& rhs) const ;
    /** Operator != */
    bool operator!=(const IndexCombination2& rhs) const ;

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& output, const IndexCombination2& out);
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
                     ParticleIndex Index3, ParticleIndex Index4);

    /** Operator < - comparison method for IndexCombination */
    bool operator< (const IndexCombination4& rhs) const ;
    /** Operator == */
    bool operator==(const IndexCombination4& rhs) const ;
    /** Operator != */
    bool operator!=(const IndexCombination4& rhs) const ;

    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& output, const IndexCombination4& out);
};

} // namespace Pomerol

#endif // #ifndef __INCLUDE_INDEX_H
