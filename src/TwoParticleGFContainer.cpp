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


#include "TwoParticleGFContainer.h"

namespace Pomerol{

TwoParticleGFContainer::TwoParticleGFContainer(const IndexClassification& IndexInfo, const StatesClassification &S,
                                               const Hamiltonian &H, const DensityMatrix &DM, const FieldOperatorContainer& Operators) :
    IndexContainer4<TwoParticleGF,TwoParticleGFContainer>(this,IndexInfo), Thermal(DM),
    S(S),H(H),DM(DM), Operators(Operators)
{}

void TwoParticleGFContainer::prepareAll(const std::set<IndexCombination4>& InitialIndices)
{
    fill(InitialIndices);
    for(std::map<IndexCombination4,ElementWithPermFreq<TwoParticleGF> >::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        static_cast<TwoParticleGF&>(iter->second).prepare();
}

void TwoParticleGFContainer::computeAll(long NumberOfMatsubaras)
{
    for(std::map<IndexCombination4,ElementWithPermFreq<TwoParticleGF> >::iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++)
        static_cast<TwoParticleGF&>(iter->second).compute(NumberOfMatsubaras);
}

TwoParticleGF* TwoParticleGFContainer::createElement(const IndexCombination4& Indices) const
{
    const AnnihilationOperator &C1 = Operators.getAnnihilationOperator(Indices.Index1);
    const AnnihilationOperator &C2 = Operators.getAnnihilationOperator(Indices.Index2);
    const CreationOperator     &CX3 = Operators.getCreationOperator   (Indices.Index3);
    const CreationOperator     &CX4 = Operators.getCreationOperator   (Indices.Index4);

    return new TwoParticleGF(S,H,C1,C2,CX3,CX4,DM);
}

std::set<IndexCombination4> TwoParticleGFContainer::getNonVanishingElements(void) const
{
    std::set<IndexCombination4> NonVanishingElements;

    for(std::map<IndexCombination4,ElementWithPermFreq<TwoParticleGF> >::const_iterator iter = ElementsMap.begin();
        iter != ElementsMap.end(); iter++){
        if(!static_cast<const TwoParticleGF&>(iter->second).isVanishing())
            NonVanishingElements.insert(iter->first);
    }
    return NonVanishingElements;
}

} // end of namespace Pomerol
