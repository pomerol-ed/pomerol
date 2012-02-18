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