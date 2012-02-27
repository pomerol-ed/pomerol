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


#include "IndexClassification.h"
#include <boost/functional/hash.hpp>

namespace Pomerol{

//
//IndexClassification::IndexInfo 
//

IndexClassification::IndexInfo::IndexInfo( const std::string &SiteLabel, const unsigned short Orbital, const unsigned short Spin):
    SiteLabel(SiteLabel), Orbital(Orbital), Spin(Spin)
{
    boost::hash<std::string> string_hash;
    SiteLabelHash = string_hash(SiteLabel);
};

bool IndexClassification::IndexInfo::operator<(const IndexClassification::IndexInfo& rhs) const 
{
    if (SiteLabelHash != rhs.SiteLabelHash) return (SiteLabelHash < rhs.SiteLabelHash);
    if (Orbital != rhs.Orbital) return (Orbital < rhs.Orbital);
    return (Spin < rhs.Spin);
}

//
//IndexClassification
//

IndexClassification::IndexClassification ( Lattice *L ) : L(L)
{
};

} // end of namespace Pomerol
