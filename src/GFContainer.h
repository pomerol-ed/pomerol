//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/


#ifndef __INCLUDE_GFCONTAINER_H
#define __INCLUDE_GFCONTAINER_H

#include"Misc.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"
#include"IndexContainer2.h"

namespace Pomerol{

typedef boost::shared_ptr<GreensFunction> GFPointer;

class GFContainer: public IndexContainer2<GreensFunction,GFContainer>, public Thermal
{

public:
    GFContainer(const IndexClassification& IndexInfo,
                const StatesClassification &S,
                const Hamiltonian &H, const DensityMatrix &DM, const FieldOperatorContainer& Operators);


    void prepareAll(const std::set<IndexCombination2>& InitialIndices = std::set<IndexCombination2>());
    void computeAll(long NumberOfMatsubaras = 0);

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
