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
#include"FourIndexObject.h"
#include"GFContainer.h"
#include"TwoParticleGFContainer.h"
#include"MatsubaraContainers.h"

namespace Pomerol{

class Vertex4Element : public Thermal {
    const TwoParticleGFContainer &Chi;
    const GFContainer &g;
    const FourIndexObject::IndexCombination& Indices;

    bool Computed;

    mutable MatsubaraContainer4<Vertex4Element> Storage;
    friend class MatsubaraContainer4<Vertex4Element>;

    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

public:

    Vertex4Element(const TwoParticleGFContainer &Chi, const GFContainer &g,
                   const FourIndexObject::IndexCombination& Indices);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    bool isComputed(void) const;
    bool isVanishing(void) const;
};

/** Objects of this class just transforms a two-particle Green's function into
 * an irreducible vertex part or into an amputated irreducible vertex.
 */
class Vertex4 : public FourIndexContainerObject<Vertex4Element>, public Thermal {

    /** A reference to a two-particle Green's function. */
    const TwoParticleGFContainer &Chi;
    /** A reference to a Green's function container */
    const GFContainer &g;

public:

    Vertex4(const IndexClassification &IndexInfo, const TwoParticleGFContainer &Chi, const GFContainer &g);

    void prepare(void);
    void prepare(const std::set<IndexCombination>& InitialCombinations);
    void compute(long NumberOfMatsubaras = 0);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_VERTEX4_H
