#include"pomerol/Index.h"

namespace Pomerol{
//
//DynamicIndexCombination
//

DynamicIndexCombination::DynamicIndexCombination(ParticleIndex N):N(N)
{
    Indices.resize(N);
    for (std::vector<ParticleIndex>::iterator it=Indices.begin(); it!=Indices.end(); ++it) *it=0;
}

DynamicIndexCombination::DynamicIndexCombination(const std::vector<ParticleIndex>& in):N(in.size()), Indices(in)
{
}

ParticleIndex& DynamicIndexCombination::operator[](const ParticleIndex position)
{
    return Indices[position];
}

const ParticleIndex DynamicIndexCombination::getIndex(const ParticleIndex position) const
{
    if (position>=N) throw (exWrongIndices());
    return Indices[position];
}

const ParticleIndex DynamicIndexCombination::getNumberOfIndices() const 
{
    return N;
}

bool DynamicIndexCombination::operator<(const DynamicIndexCombination& rhs) const
{
    if (rhs.N!=N) throw (exWrongIndices());
    bool result = Indices[0] < rhs.Indices[0];
    for (ParticleIndex i=0; i<N-1 && Indices[i]==rhs.Indices[i]; ++i) result=(Indices[i+1]<rhs.Indices[i+1]);
    return result;
}

bool DynamicIndexCombination::operator==(const DynamicIndexCombination& rhs) const 
{
    bool result = true;
    for (ParticleIndex i=0; i<N && result; ++i) result=(Indices[i]==rhs.Indices[i]);
    return result;
}

bool DynamicIndexCombination::operator!=(const DynamicIndexCombination& rhs) const 
{
    return !((*this)==rhs);
}

DynamicIndexCombination& DynamicIndexCombination::operator=(const DynamicIndexCombination& rhs) 
{
    if (this != &rhs) {
        N=rhs.N;
        Indices=rhs.Indices;
    };
    return (*this);
}

std::ostream& operator<<(std::ostream& output, const DynamicIndexCombination& out)
{
    output << "("; 
    for (std::vector<ParticleIndex>::const_iterator it1=out.Indices.begin(); it1!=out.Indices.end(); ++it1) 
        output << *it1; 
    output << ")" << std::flush;
    return output;
}


const char* DynamicIndexCombination::exWrongIndices::what() const throw(){
    return "Wrong indices";
};
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
