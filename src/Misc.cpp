//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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

/** \file Misc.cpp
**  \brief Various useful code.
** 
** \author    Igor Krivenko (igor@shg.ru)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#include"Misc.h"

namespace Pomerol {
    //////////////////
// Permutation3 //
//////////////////
bool Permutation3::operator==(const Permutation3& rhs) const
{
    return (sign==rhs.sign && perm[0] == rhs.perm[0] && perm[1]==rhs.perm[1]);
}

bool Permutation3::operator!=(const Permutation3& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& out, const Permutation3 &rhs)
{
    out << (rhs.sign==-1?"-":" ") << rhs.perm[0]+1 << rhs.perm[1]+1 << rhs.perm[2]+1 << std::flush;
    return out;
}

const Permutation3 permutations3[6] = {
    {{0,1,2},1},
    {{0,2,1},-1},
    {{1,0,2},-1},
    {{1,2,0},1},
    {{2,0,1},1},
    {{2,1,0},-1}
};

//////////////////
// Permutation4 //
//////////////////
bool Permutation4::operator==(const Permutation4& rhs) const
{
    return (sign==rhs.sign && perm[0] == rhs.perm[0] && perm[1]==rhs.perm[1] && perm[2] == rhs.perm[2]);
}

bool Permutation4::operator!=(const Permutation4& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& out, const Permutation4 &p)
{
    out << (p.sign==-1?"-":" ") << p.perm[0]+1 << p.perm[1]+1 << p.perm[2]+1 << p.perm[3]+1 << std::flush;
    return out;
}

const Permutation4 permutations4[24] = {
    {{0,1,2,3}, 1},  {{0,1,3,2},-1},  {{0,2,1,3},-1},  {{0,2,3,1}, 1},  {{0,3,1,2}, 1},  {{0,3,2,1},-1},
    {{1,0,2,3},-1},  {{1,0,3,2}, 1},  {{1,2,0,3}, 1},  {{1,2,3,0},-1},  {{1,3,0,2},-1},  {{1,3,2,0}, 1},
    {{2,0,1,3}, 1},  {{2,0,3,1},-1},  {{2,1,0,3},-1},  {{2,1,3,0}, 1},  {{2,3,0,1}, 1},  {{2,3,1,0},-1},
    {{3,0,1,2},-1},  {{3,0,2,1}, 1},  {{3,1,0,2}, 1},  {{3,1,2,0},-1},  {{3,2,0,1},-1},  {{3,2,1,0}, 1}
};

} 