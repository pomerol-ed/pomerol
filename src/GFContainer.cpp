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


#include "GFContainer.h"

namespace Pomerol{

GFContainer::IndexCombination::IndexCombination(ParticleIndex Index1, ParticleIndex Index2) :
    Index1(Index1), // Index of C
    Index2(Index2) // Index of C^+
{}

bool GFContainer::IndexCombination::operator<(const GFContainer::IndexCombination& rhs) const
{
    return (Index1<rhs.Index1) || (Index1==rhs.Index1 && Index2 < rhs.Index2);
}


std::ostream& operator<<(std::ostream& output,const GFContainer::IndexCombination& out)
{
    output << "(" << out.Index1 << out.Index2 << ")" << std::flush;
    return output;
}

/*=========================================================================*/

GFContainer::GFContainer ( StatesClassification &S, const Hamiltonian &H, const DensityMatrix &DM,
                           IndexClassification& IndexInfo, const FieldOperatorContainer& Operators) :
    Thermal(DM),S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators)
{}

void GFContainer::prepare(void)
{
    // If no initial index combinations are supplied, enumerate all possible combinations
    for (ParticleIndex Index1=0; Index1<IndexInfo.getIndexSize(); ++Index1)
    for (ParticleIndex Index2=0; Index2<IndexInfo.getIndexSize(); ++Index2){
        GreensFunction *pGF = new GreensFunction (S,H,Operators.getAnnihilationOperator(Index1),Operators.getCreationOperator(Index2),DM);
        pGF->prepare();
        if(!pGF->isVanishing())
            mapGreensFunctions[IndexCombination(Index1,Index2)] = pGF;
    }
}

void GFContainer::prepare(const std::vector<IndexCombination*>& InitialIndices)
{
    for (std::vector<IndexCombination*>::const_iterator it=InitialIndices.begin(); it!=InitialIndices.end(); ++it) {
        GreensFunction *GF = new GreensFunction (S,H,Operators.getAnnihilationOperator((*it)->Index1),Operators.getCreationOperator((*it)->Index2),DM);
        GF->prepare();
        if (!GF->isVanishing())
            mapGreensFunctions[**it]=GF;
    }
}

void GFContainer::computeValues(long NumberOfMatsubaras)
{
    for (std::map<IndexCombination,GreensFunction*>::iterator it=mapGreensFunctions.begin();it!=mapGreensFunctions.end();++it){
        INFO("GFContainer: computing G_{" << it->first << "}");
        it->second->precomputeParts();
        it->second->computeValues(NumberOfMatsubaras);
    }
}

const MatrixType& GFContainer::operator()(long MatsubaraNumber) const
{
    MatrixType* Output = new MatrixType(IndexInfo.getIndexSize(), IndexInfo.getIndexSize());
    Output->setZero();
    for (std::map<IndexCombination, GreensFunction*>::const_iterator it1=mapGreensFunctions.begin(); it1!=mapGreensFunctions.end(); it1++){
        (*Output)(it1->first.Index1,it1->first.Index2) = (*it1->second)(MatsubaraNumber);
    }
    return *Output; 
}

ComplexType GFContainer::operator()(ParticleIndex Index1, ParticleIndex Index2, long MatsubaraNumber) const
{
    return (!isVanishing(Index1,Index2)) ? 
        (*(mapGreensFunctions.find(IndexCombination(Index1,Index2))->second))(MatsubaraNumber) : 
        ((Index1==Index2) ? 1.0/(static_cast<RealType>(2*MatsubaraNumber+1)*MatsubaraSpacing) : 0.);
}

bool GFContainer::isVanishing(ParticleIndex Index1, ParticleIndex Index2) const
{
   return (mapGreensFunctions.count(IndexCombination(Index1,Index2))==0); 
}

} // end of namespace Pomerol
