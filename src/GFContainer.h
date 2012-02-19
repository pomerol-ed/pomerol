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


/** \file src/GFContainer.h
**
** \author Andrey Antipov (antipov@ct-qmc.org)
** \author Igor Krivenko (igor@shg.ru)
*/


#ifndef __INCLUDE_GFCONTAINER_H
#define __INCLUDE_GFCONTAINER_H

#include"Misc.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"
#include"IndexContainer2.h"

#include<set>

namespace Pomerol{

class GFContainer: public IndexContainer2<GreensFunction,GFContainer>, public Thermal
{
    typedef boost::shared_ptr<GreensFunction> GFPointer;
    typedef std::pair<IndexCombination2,GFPointer> IndicesGFPair;

public:
    GFContainer(const IndexClassification& IndexInfo,
                const StatesClassification &S,
                const Hamiltonian &H, const DensityMatrix &DM, const FieldOperatorContainer& Operators);


    void prepareAll(const std::set<IndexCombination2>& InitialIndices = std::set<IndexCombination2>());
    void computeAll(long NumberOfMatsubaras = 0);

    std::set<IndexCombination2> getNonVanishingElements(void) const;

protected:

    friend class IndexContainer2<GreensFunction,GFContainer>;
    GreensFunction* createElement(const IndexCombination2& Indices) const;

    const StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM;
    const FieldOperatorContainer &Operators;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
