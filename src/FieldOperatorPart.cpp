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


#include <fstream>

#include "FieldOperatorPart.h"

using std::stringstream;

namespace Pomerol{

FieldOperatorPart::FieldOperatorPart(
        IndexClassification &IndexInfo, StatesClassification &S, HamiltonianPart &h_from,  HamiltonianPart &h_to, ParticleIndex i) : 
        ComputableObject(), IndexInfo(IndexInfo), S(S), h_from(h_from), h_to(h_to), i(i)
{};


ColMajorMatrixType& FieldOperatorPart::getColMajorValue()
{
    return elementsColMajor;
}

RowMajorMatrixType& FieldOperatorPart::getRowMajorValue()
{
    return elementsRowMajor;
}

void FieldOperatorPart::print_to_screen()  //print to screen C and CX
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    for (int P=0; P<elementsColMajor.outerSize(); ++P)
        for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it)
        {
                QuantumState N = S.getQuantumStates(to)[it.row()];
                QuantumState M = S.getQuantumStates(from)[it.col()];
                std::cout << N <<" " << M << " : " << it.value() << std::endl;
        };
}

void FieldOperatorPart::dump() //writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
DEBUG("No dump method is implemented for field operators");
assert(0);
/*    std::stringstream filename;
    filename << (*this).OUT.fullpath() << "/" << "part" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
    ofstream outCpart;
    outCpart.open(filename.str().c_str());
     for (int P=0; P<elementsColMajor.outerSize(); ++P)
        for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it)
        {
                QuantumState N = it.row();//S.getQuantumStates(to)[it.row()];
                QuantumState M = it.col();//S.getQuantumStates(from)[it.col()];
                outCpart << S.getQuantumStates(h_to.id())[N] <<" " << S.getQuantumStates(h_from.id())[M] << "  " << it.value() << endl;
        };


 //   outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements.toDense() << endl;

    outCpart.close();
    cout << "The part of field operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
*/
}

void FieldOperatorPart::compute()
{
    QuantumNumbers to = h_to.id();
    QuantumNumbers from = h_from.id();

    const std::vector<QuantumState>& toStates = S.getQuantumStates(to);
    const std::vector<QuantumState>& fromStates = S.getQuantumStates(from);
    
    DynamicSparseMatrixType tempElements(toStates.size(),fromStates.size());
    
    for (std::vector<QuantumState>::const_iterator current_state = toStates.begin();
                                                   current_state < toStates.end(); current_state++)
    {     
        QuantumState L=*current_state;
        if (checkL(L))
        {
            int K = retK(L);
            if( (mFunc(L,K,i)!= 0) )
            {                                                 
                int l=S.getInnerState(L), k=S.getInnerState(K);                             // l,k in part of Hamilt                        
               
                for ( unsigned int n=0; n<toStates.size(); n++)
                {
                    if(h_to.reH(l,n)!=0)
                    {
                        for (unsigned int m=0; m<fromStates.size(); m++)
                        {
                            RealType C_nm = h_to.reH(l,n)*mFunc(L,K,i)*h_from.reH(k,m);
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

// Functions of specialized classes

AnnihilationOperatorPart::AnnihilationOperatorPart(IndexClassification &IndexInfo, StatesClassification &S, 
                                                   HamiltonianPart &h_from, HamiltonianPart &h_to, ParticleIndex i) :
FieldOperatorPart(IndexInfo,S,h_from,h_to,i)
{}

CreationOperatorPart::CreationOperatorPart(IndexClassification &IndexInfo, StatesClassification &S, 
                                                  HamiltonianPart &h_from, HamiltonianPart &h_to, ParticleIndex i) :
FieldOperatorPart(IndexInfo,S,h_from,h_to,i)
{}

QuantumState AnnihilationOperatorPart::retK(QuantumState L)							//return K for C
{	
	return( L + (1<<i) );
}

QuantumState CreationOperatorPart::retK(QuantumState L)							//return K for CX
{	

	return( L - (1<<i) );

}

int AnnihilationOperatorPart::mFunc(QuantumState state1, QuantumState state2, ParticleIndex i)
{
	int flag=1, p=0;
	for (int m=0; m<IndexInfo.getIndexSize(); m++)
	{
		if( m == i )
		{
			if ((S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)
	{
		for (int m=0;m<i;m++) p+=S.n_i(state2,m);
	}
	return (flag*(1-2*(p%2)));
}

int CreationOperatorPart::mFunc(QuantumState state1, QuantumState state2, ParticleIndex i)
{
	int flag=1, p=0;
	for (int m=0; m<IndexInfo.getIndexSize(); m++)
	{
		if( m == i )
		{
			if ((S.n_i(state2,i)==1) || (S.n_i(state1,i)==0) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)
	{
		for (int m=0;m<i;m++) p+=S.n_i(state2,m);
	}
	return (flag*(1-2*(p%2)));
}


bool AnnihilationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==0);
}

bool CreationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==1);
}

BlockNumber FieldOperatorPart::getLeftIndex()
{
	return h_to.getId();
};

BlockNumber FieldOperatorPart::getRightIndex()
{
	return h_from.getId();
};

CreationOperatorPart& AnnihilationOperatorPart::transpose()
{
	CreationOperatorPart *CX = new CreationOperatorPart(IndexInfo, S, h_to, h_from, i); // swapped h_to and h_from
	CX->elementsRowMajor = elementsRowMajor.transpose();
	CX->elementsColMajor = elementsColMajor.transpose();
	return *CX;
};

AnnihilationOperatorPart& CreationOperatorPart::transpose()
{
	AnnihilationOperatorPart *C = new AnnihilationOperatorPart(IndexInfo, S, h_to, h_from, i); // swapped h_to and h_from
	C->elementsRowMajor = elementsRowMajor.transpose();
	C->elementsColMajor = elementsColMajor.transpose();
	return *C;
};

} // end of namespace Pomerol
