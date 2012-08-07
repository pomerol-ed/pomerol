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
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>


#include"Index.h"

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
