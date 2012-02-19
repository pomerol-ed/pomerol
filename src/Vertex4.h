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


/** \file src/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_VERTEX4_H
#define __INCLUDE_VERTEX4_H

#include"Misc.h"
#include"Logger.h"
#include"GFContainer.h"
#include"TwoParticleGFContainer.h"
#include"MatsubaraContainers.h"

namespace Pomerol{

class Vertex4 : public Thermal, public ComputableObject {

    enum {Constructed,Computed};

    TwoParticleGF &Chi4;
    GreensFunction &G13;
    GreensFunction &G24;
    GreensFunction &G14;
    GreensFunction &G23;

    /** Storage for precomputed values. */
    mutable MatsubaraContainer4<Vertex4> Storage;
    friend class MatsubaraContainer4<Vertex4>;

    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

public:

    Vertex4(TwoParticleGF& Chi4,
            GreensFunction& G13, GreensFunction& G24,
            GreensFunction& G14, GreensFunction& G23);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    bool isVanishing(void) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_VERTEX4_H
