#include "pomerol/Misc.h"

namespace Pomerol {

//////////
// spin //
//////////
std::ostream & operator<<(std::ostream & os, spin s)
{
    return os << (s == up ? "up" : "dn");
}


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

} // namespace Pomerol
