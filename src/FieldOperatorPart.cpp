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


#include "FieldOperatorPart.h"

using std::stringstream;

namespace Pomerol{

FieldOperatorPart::FieldOperatorPart(
        const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom,  const HamiltonianPart &HTo, ParticleIndex PIndex) : 
        IndexInfo(IndexInfo), S(S), HFrom(HFrom), HTo(HTo), PIndex(PIndex)
{}

void FieldOperatorPart::compute()
{
    QuantumNumbers to = HTo.getQuantumNumbers();
    QuantumNumbers from = HFrom.getQuantumNumbers();
 
    const std::vector<QuantumState>& toStates = S.getQuantumStates(to);
    const std::vector<QuantumState>& fromStates = S.getQuantumStates(from);

    DynamicSparseMatrixType tempElements(toStates.size(),fromStates.size());

    for (std::vector<QuantumState>::const_iterator CurrentState = toStates.begin();
                                                   CurrentState < toStates.end(); CurrentState++)
    {
	QuantumState L=*CurrentState;
	if (checkL(L))
	{
	    int K = retK(L);
	    if( (mFunc(L,K,PIndex)!= 0) )
	    {
		int l=S.getInnerState(L), k=S.getInnerState(K);                             // l,k in part of Hamilt                        
		for (size_t n=0; n<toStates.size(); n++)
		{
		    if(HTo.getMatrixElement(l,n)!=0)
		    {
			for (size_t m=0; m<fromStates.size(); m++)
			{
			    RealType C_nm = HTo.getMatrixElement(l,n)*mFunc(L,K,PIndex)*HFrom.getMatrixElement(k,m);
			    if (fabs(C_nm)>MatrixElementTolerance)
			    {
				tempElements.coeffRef(n,m) += C_nm;
			    }
			}
		    }
		}
	    }
	}
    }
    tempElements.prune(MatrixElementTolerance);
    elementsRowMajor = tempElements;
    elementsColMajor = tempElements;
}

const ColMajorMatrixType& FieldOperatorPart::getColMajorValue(void) const
{
    return elementsColMajor;
}
 
const RowMajorMatrixType& FieldOperatorPart::getRowMajorValue(void) const
{
    return elementsRowMajor;
}

void FieldOperatorPart::print_to_screen()  //print to screen C and CX
{
    QuantumNumbers to   = HTo.getQuantumNumbers();
    QuantumNumbers from = HFrom.getQuantumNumbers();
    for (size_t P=0; P<elementsColMajor.outerSize(); ++P)
	for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it)
        {
	    QuantumState N = S.getQuantumStates(to)[it.row()];
	    QuantumState M = S.getQuantumStates(from)[it.col()];
	    std::cout << N <<" " << M << " : " << it.value() << std::endl;
        }
}

BlockNumber FieldOperatorPart::getLeftIndex(void) const
{
    return HTo.getBlockNumber();
}

BlockNumber FieldOperatorPart::getRightIndex(void) const
{
    return HFrom.getBlockNumber();
}

// Functions of specialized classes

AnnihilationOperatorPart::AnnihilationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, 
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{}

CreationOperatorPart::CreationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, 
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{}

QuantumState AnnihilationOperatorPart::retK(QuantumState L) const //return K for C
{
    return(L + (1<<PIndex));
}

QuantumState CreationOperatorPart::retK(QuantumState L) const	//return K for CX
{
    return( L - (1<<PIndex) );
}

int AnnihilationOperatorPart::mFunc(QuantumState state1, QuantumState state2, ParticleIndex PIndex) const
{
    int flag=1, p=0;
    for (int m=0; m<IndexInfo.getIndexSize(); m++)
    {
	if( m == PIndex )
	{
	    if ((S.n_i(state2,PIndex)==0) || (S.n_i(state1,PIndex)==1) ) {flag=0;break;} else flag=1;
	}
	else
	{
	    if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
	}
    }
    if (flag==1)
    {
	for (int m=0;m<PIndex;m++) p+=S.n_i(state2,m);
    }
    return (flag*(1-2*(p%2)));
}

int CreationOperatorPart::mFunc(QuantumState state1, QuantumState state2, ParticleIndex PIndex) const
{
    int flag=1, p=0;
    for (int m=0; m<IndexInfo.getIndexSize(); m++)
    {
	if( m == PIndex )
	{
	    if ((S.n_i(state2,PIndex)==1) || (S.n_i(state1,PIndex)==0) ) {flag=0;break;} else flag=1;
	}
	else
	{
	    if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
	}
    }
    if (flag==1)
    {
	for (int m=0;m<PIndex;m++) p+=S.n_i(state2,m);
    }
    return (flag*(1-2*(p%2)));
}

bool AnnihilationOperatorPart::checkL(QuantumState L) const
{
    return (S.n_i(L,PIndex)==0);
}

bool CreationOperatorPart::checkL(QuantumState L) const
{
    return (S.n_i(L,PIndex)==1);
}

const CreationOperatorPart& AnnihilationOperatorPart::transpose() const
{
    CreationOperatorPart *CX = new CreationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    CX->elementsRowMajor = elementsRowMajor.transpose();
    CX->elementsColMajor = elementsColMajor.transpose();
    return *CX;
}

const AnnihilationOperatorPart& CreationOperatorPart::transpose() const
{
    AnnihilationOperatorPart *C = new AnnihilationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    C->elementsRowMajor = elementsRowMajor.transpose();
    C->elementsColMajor = elementsColMajor.transpose();
    return *C;
}

} // end of namespace Pomerol
