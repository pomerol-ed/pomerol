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


#ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
#define __INCLUDE_TWOPARTICLEGFCONTAINER_H

#include"Misc.h"
#include"TwoParticleGF.h"
#include"FieldOperatorContainer.h"
#include"IndexContainer4.h"

namespace Pomerol{

typedef boost::shared_ptr<TwoParticleGF> GF2Pointer;

class TwoParticleGFContainer: public IndexContainer4<TwoParticleGF,TwoParticleGFContainer>, public Thermal
{
public:

    TwoParticleGFContainer(const IndexClassification& IndexInfo, const StatesClassification &S,
                           const Hamiltonian &H, const DensityMatrix &DM, const FieldOperatorContainer& Operators);

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
#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGFCONTAINER_H
