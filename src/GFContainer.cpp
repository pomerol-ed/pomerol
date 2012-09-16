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


#include "GFContainer.h"

namespace Pomerol{

GFContainer::GFContainer ( const IndexClassification& IndexInfo,
                           const StatesClassification& S,
                           const Hamiltonian &H, const DensityMatrix &DM,
                           const FieldOperatorContainer& Operators) :
    IndexContainer2<GreensFunction,GFContainer>(this,IndexInfo),
    Thermal(DM), S(S), H(H), DM(DM), Operators(Operators)
{}

void GFContainer::prepareAll(const std::set<IndexCombination2>& InitialIndices)
{
    fill(InitialIndices);
    for(std::map<IndexCombination2,GFPointer>::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        (iter->second)->prepare();
}

void GFContainer::computeAll(long NumberOfMatsubaras)
{
    for(std::map<IndexCombination2,GFPointer>::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        (iter->second)->compute(NumberOfMatsubaras);
}

GreensFunction* GFContainer::createElement(const IndexCombination2& Indices) const
{
    return new GreensFunction(S,H, Operators.getAnnihilationOperator(Indices.Index1),
                                   Operators.getCreationOperator(Indices.Index2),DM);
}

} // end of namespace Pomerol
