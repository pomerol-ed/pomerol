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


#include "GFContainer.h"
extern std::ostream& OUTPUT_STREAM;

namespace Pomerol{

GFContainer::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cdagindex2)
{
    Indices[0]=cindex1;
    Indices[1]=cdagindex2;
}

bool GFContainer::IndexCombination::operator<(const GFContainer::IndexCombination& rhs) const
{
  return (Indices[0]<rhs.Indices[0]) || (Indices[0]==rhs.Indices[0] && Indices[1] < rhs.Indices[1]);
}


std::ostream& operator<<(std::ostream& output,const GFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << ")" << std::flush;
return output;
}


/*=========================================================================*/

GFContainer::GFContainer ( StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,IndexClassification& IndexInfo, FieldOperatorContainer& Operators):
  ComputableObject(),Thermal(DM),S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators)
{
};

MatrixType& GFContainer::operator()(long MatsubaraNumber)
{
    MatrixType* Output = new MatrixType(IndexInfo.getIndexSize(), IndexInfo.getIndexSize());
    Output->setZero();
    for (std::map<IndexCombination, GreensFunction*>::iterator it1=mapGreensFunctions.begin(); it1!=mapGreensFunctions.end(); it1++){
        (*Output)(it1->first.Indices[0],it1->first.Indices[1]) = (*it1->second)(MatsubaraNumber);
        };
    return *Output; 
};

ComplexType GFContainer::operator()(ParticleIndex i, ParticleIndex j, long MatsubaraNumber)
{
    return (!vanishes(i,j))?(*mapGreensFunctions[IndexCombination(i,j)])(MatsubaraNumber):((i==j)?I*(-1.)*beta/((2*MatsubaraNumber+1)*M_PI):0.);
    //return (!vanishes(i,j))?(*mapGreensFunctions[IndexCombination(i,j)])(MatsubaraNumber):((i==j)?1.:0.);
};

bool GFContainer::vanishes(ParticleIndex i, ParticleIndex j)
{
   return (mapGreensFunctions.count(IndexCombination(i,j))==0); 
};

void GFContainer::defineInitialIndices()
{
    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i)
        for (ParticleIndex j=0; j<IndexInfo.getIndexSize(); ++j){
            IndexCombination * Gij = new IndexCombination(i,j);
            InitialIndices.push_back(Gij);
            };
}

void GFContainer::readInitialIndices(std::vector<IndexCombination*> &in)
{
    InitialIndices=in;
};

void GFContainer::prepare()
{
if (Status == Constructed){
    if (InitialIndices.size()==0) defineInitialIndices();
    for (std::vector<IndexCombination*>::const_iterator it=InitialIndices.begin(); it!=InitialIndices.end(); ++it) {
        GreensFunction *GF = new GreensFunction (S,H,Operators.getAnnihilationOperator((*it)->Indices[0]),Operators.getCreationOperator((*it)->Indices[1]),DM);
        GF->prepare();
        if (!GF->vanishes()) mapGreensFunctions[**it]=GF;
        }
    Status = Prepared;
    };
};

void GFContainer::compute()
{
if (Status == Prepared){
    for (std::map<IndexCombination,GreensFunction*>::iterator it1=mapGreensFunctions.begin();it1!=mapGreensFunctions.end();++it1){
           INFO("GFContainer: computing G_{" << it1->first << "}");
           it1->second->compute(); 
        };
    Status = Computed;
    };
};

void GFContainer::dumpToPlainText(long wn)
{
    for (std::map<IndexCombination,GreensFunction*>::iterator it1=mapGreensFunctions.begin();it1!=mapGreensFunctions.end();++it1){
           it1->second->dumpToPlainText(wn); 
        };
};

} // end of namespace Pomerol
