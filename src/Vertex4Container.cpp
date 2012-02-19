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


#include "Vertex4Container.h"

namespace Pomerol{

Vertex4Container::Vertex4Container(TwoParticleGFContainer& Chi4, GFContainer& G) :
    IndexContainer4(this,IndexInfo), Thermal(Chi4),
    Chi4(Chi4), G(G)
{}

void Vertex4Container::prepareAll(const std::set<IndexCombination4>& InitialIndices)
{
    fill(InitialIndices);
}

void Vertex4Container::computeAll(long NumberOfMatsubaras)
{
    for(std::map<IndexCombination4,ElementWithPermFreq<Vertex4> >::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        static_cast<Vertex4&>(iter->second).compute(NumberOfMatsubaras);
}

Vertex4* Vertex4Container::createElement(const IndexCombination4& Indices) const
{
    GreensFunction& G13 = G(Indices.Index1,Indices.Index3);
    GreensFunction& G24 = G(Indices.Index2,Indices.Index4);
    GreensFunction& G14 = G(Indices.Index1,Indices.Index4);
    GreensFunction& G23 = G(Indices.Index2,Indices.Index3);

    return new Vertex4(Chi4(Indices),G13,G24,G14,G23);
}

} // end of namespace Pomerol
