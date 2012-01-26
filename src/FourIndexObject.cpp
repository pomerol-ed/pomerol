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
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#include "FourIndexObject.h"
//
//FourIndexObject::IndexCombination
//

namespace Pomerol{

FourIndexObject::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4)
{
    Indices[0]=cindex1;
    Indices[1]=cindex2;
    Indices[2]=cdagindex3;
    Indices[3]=cdagindex4;
}

bool FourIndexObject::IndexCombination::operator<(const FourIndexObject::IndexCombination& rhs) const
{
  return (Indices[0] < rhs.Indices[0]) || 
         (Indices[0] == rhs.Indices[0] && Indices[1] < rhs.Indices[1] ) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] < rhs.Indices[2]) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] < rhs.Indices[3]); 
}

bool FourIndexObject::IndexCombination::operator==(const FourIndexObject::IndexCombination& rhs) const
{
    return (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] == rhs.Indices[3]); 
}

bool FourIndexObject::IndexCombination::operator!=(const FourIndexObject::IndexCombination& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& output,const FourIndexObject::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

//
// FourIndexContainerObject
//
const Permutation4 FourIndexContainerObject::TrivialOperatorPermutations[4] = { permutations4[0], permutations4[1], permutations4[6], permutations4[7] };

} // end of namespace Pomerol
