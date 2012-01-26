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
** \brief A container for either creation or annihilation operators in eigenvector basis
**
** \author Andrey Antipov (antipov@ct-qmc.org)
*/


#ifndef __INCLUDE_GFCONTAINER_H
#define __INCLUDE_GFCONTAINER_H

#include"Misc.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"

namespace Pomerol{

class GFContainer : public Thermal 
{

public:
    struct IndexCombination; 
    GFContainer(StatesClassification &S, const Hamiltonian &H, const DensityMatrix &DM, IndexClassification& IndexInfo, const FieldOperatorContainer& Operators);

    void prepare(void);
    void prepare(const std::vector<IndexCombination*>& InitialIndices);
    void computeValues(long NumberOfMatsubaras);

    bool isVanishing(ParticleIndex Index1, ParticleIndex Index2) const;

    const MatrixType& operator()(long MatsubaraNumber) const;
    ComplexType operator()(ParticleIndex Index1, ParticleIndex Index2, long MatsubaraNumber) const;

private:
    IndexClassification &IndexInfo;
    StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM;
    const FieldOperatorContainer &Operators;

    std::map<IndexCombination, GreensFunction*> mapGreensFunctions;
};

struct GFContainer::IndexCombination
{
    const ParticleIndex Index1, Index2;

    bool operator<(const GFContainer::IndexCombination& rhs) const;
    IndexCombination(ParticleIndex Index1, ParticleIndex Index2);

    friend std::ostream& operator<<(std::ostream& output, const GFContainer::IndexCombination& out);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
