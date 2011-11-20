//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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
#include"ComputableObject.h"
#include"GreensFunction.h"
#include"FieldOperatorContainer.h"

namespace Pomerol{

class GFContainer : public ComputableObject, public Thermal 
{

public:
    struct IndexCombination; 
    GFContainer(StatesClassification &S, Hamiltonian &H, DensityMatrix &DM, IndexClassification& IndexInfo, FieldOperatorContainer& Operators);
    void readInitialIndices(std::vector<IndexCombination*> &in);
    void prepare();
    void compute();
    bool vanishes(ParticleIndex i, ParticleIndex j);
    MatrixType& operator()(long MatsubaraNumber);
    ComplexType operator()(ParticleIndex i, ParticleIndex j, long MatsubaraNumber);
    void dumpToPlainText(long wn);
    std::vector<IndexCombination*> InitialIndices;
private:
    StatesClassification &S;
    Hamiltonian &H;
    DensityMatrix &DM;
    IndexClassification &IndexInfo;
    FieldOperatorContainer &Operators;

    std::map<IndexCombination, GreensFunction*> mapGreensFunctions;
    void defineInitialIndices();
};

struct GFContainer::IndexCombination
{
    ParticleIndex Indices[2];
    friend std::ostream& operator<<(std::ostream& output, const GFContainer::IndexCombination& out);
    bool operator<(const GFContainer::IndexCombination& rhs) const ;
    IndexCombination(ParticleIndex cindex1, ParticleIndex cdagindex2);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_GFCONTAINER_H
