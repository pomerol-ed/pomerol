/** \file Index.h
**  \brief Types related to indices and their combinations.
** 
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef __INCLUDE_INDEX_H
#define __INCLUDE_INDEX_H

#include"Misc.h"

namespace Pomerol {

template <int N> struct StaticIndexCombination
{
/** A structure to handle a combination of 2 Particle Indices. */
    const ParticleIndex Indices[N];

    /** Constructor 
     * \param[in] in 1d array with N ParticleIndex numbers.
     */
    StaticIndexCombination (ParticleIndex in [][1][N]);
    /** Operator < - comparison method for IndexCombination */
    bool operator < (const StaticIndexCombination<N>& rhs) const ;
    /** Operator == */
    bool operator==(const StaticIndexCombination<N>& rhs) const ;
    /** Operator != */
    bool operator!=(const StaticIndexCombination<N>& rhs) const ;

    /** Output to external stream */
    template <int M>
    friend std::ostream& operator<< (std::ostream& output, const StaticIndexCombination<M>& out);

};

/** A structure to hold a combination of several indices. This is an implementation for a dynamic number of Indices. Slower than a static one. */
struct DynamicIndexCombination 
{
protected:
    /** Total number of ParticleIndices. */
    ParticleIndex N;
    /** Indices. */
    std::vector<ParticleIndex> Indices;
public:
    /** Operator < - comparison method for IndexCombination */
    bool operator < (const DynamicIndexCombination& rhs) const ;
    /** Operator == */
    bool operator==(const DynamicIndexCombination& rhs) const ;
    /** Operator != */
    bool operator!=(const DynamicIndexCombination& rhs) const ;
    /** Operator = */
    DynamicIndexCombination& operator=(const DynamicIndexCombination& rhs);

    /** Output to external stream */
    friend std::ostream& operator<< (std::ostream& output, const DynamicIndexCombination& out);

    /** Returns total number of indices in current combination. */
    const ParticleIndex getNumberOfIndices() const;
    /** Returns index at given position. Safe. */
    const ParticleIndex getIndex(const ParticleIndex position) const;
    /** Returns index at given position. It can be changed. Unsafe. */
    ParticleIndex& operator[](const ParticleIndex position);

    /** Constructor. Sets all indices to 0.
     * \param[in] N Number of indices. 
     */
    DynamicIndexCombination(ParticleIndex N);

    /** Constructor. Sets all indices to provided values.
     * \param[in] in A vector of indices. 
     */
    DynamicIndexCombination(const std::vector<ParticleIndex>& in);

    class exWrongIndices : public std::exception { virtual const char* what() const throw(); };
};

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

} // end of namespace Pomerol

#endif // #ifndef __INCLUDE_INDEX_H
