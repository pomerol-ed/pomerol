#include"pomerol/Index.h"

namespace Pomerol{

///////////////////////
// IndexCombination2 //
///////////////////////

IndexCombination2::IndexCombination2(ParticleIndex Index1, ParticleIndex Index2) :
    Index1(Index1), // Index of C
    Index2(Index2) // Index of C^+
{}

bool IndexCombination2::operator<(const IndexCombination2& rhs) const
{
    return (Index1<rhs.Index1) || (Index1==rhs.Index1 && Index2 < rhs.Index2);
}

std::ostream& operator<<(std::ostream& output, const IndexCombination2& out)
{
    output << "(" << out.Index1 << out.Index2 << ")" << std::flush;
    return output;
}

///////////////////////
// IndexCombination4 //
///////////////////////

IndexCombination4::IndexCombination4(ParticleIndex Index1, ParticleIndex Index2,
                                     ParticleIndex Index3, ParticleIndex Index4) :
    Index1(Index1), Index2(Index2), Index3(Index3), Index4(Index4)
{}

bool IndexCombination4::operator<(const IndexCombination4& rhs) const
{
  return (Index1 < rhs.Index1) ||
         (Index1 == rhs.Index1 && Index2 < rhs.Index2 ) ||
         (Index1 == rhs.Index1 && Index2 == rhs.Index2 && Index3 < rhs.Index3) ||
         (Index1 == rhs.Index1 && Index2 == rhs.Index2 && Index3 == rhs.Index3 && Index4 < rhs.Index4);
}

bool IndexCombination4::operator==(const IndexCombination4& rhs) const
{
    return (Index1 == rhs.Index1 && Index2 == rhs.Index2 && Index3 == rhs.Index3 && Index4 == rhs.Index4);
}

bool IndexCombination4::operator!=(const IndexCombination4& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& output,const IndexCombination4& out)
{
    output << "(" << out.Index1 << out.Index2 << out.Index3 << out.Index4 << ")";
    return output;
}

} // end of namespace Pomerol
