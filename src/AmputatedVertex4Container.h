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


#ifndef __INCLUDE_AMPUTATEDVERTEX4CONTAINER_H
#define __INCLUDE_AMPUTATEDVERTEX4CONTAINER_H

#include"Misc.h"
#include"GFContainer.h"
#include"AmputatedVertex4.h"
#include"IndexContainer4.h"

namespace Pomerol{

typedef boost::shared_ptr<AmputatedVertex4> AV4Pointer;

class AmputatedVertex4Container : public IndexContainer4<AmputatedVertex4,AmputatedVertex4Container>, public Thermal
{
public:

    AmputatedVertex4Container(const IndexClassification& IndexInfo, const StatesClassification &S,
                           const TwoParticleGFContainer& Chi4, const GFContainer& G);

    void prepareAll(const std::set<IndexCombination4>& InitialIndices = std::set<IndexCombination4>());
    void computeAll(long NumberOfMatsubaras = 0);

protected:

    friend class IndexContainer4<TwoParticleGF,TwoParticleGFContainer>;
    TwoParticleGF* createElement(const IndexCombination4& Indices) const;

    const StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM; 
    const FieldOperatorContainer &Operators;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_AMPUTATEDVERTEX4CONTAINER_H
