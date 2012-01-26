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


#include"HamiltonianPart.h"
#include"StatesClassification.h"
#include<sstream>
#include<Eigen/Eigenvalues>

// class HamiltonianPart

namespace Pomerol{

HamiltonianPart::HamiltonianPart(IndexClassification& F, StatesClassification& S, const QuantumNumbers& QN):
    IndexInfo(F), S(S), QN(QN)
{}

void HamiltonianPart::add_nTerm(InnerQuantumState InnerState, nTerm *N)
{
    QuantumState state = S.getQuantumState(QN,InnerState);	//real state	
	H(InnerState,InnerState) += N->Value*S.n_i(state,N->bit[0]);
}

void HamiltonianPart::add_nnTerm(InnerQuantumState InnerState, nnTerm *T)
{
    QuantumState state = S.getQuantumState(QN,InnerState);	//real state
	H(InnerState,InnerState) += T->Value*S.n_i(state,T->bit[0])*S.n_i(state,T->bit[2]);
}

void HamiltonianPart::add_spinflipTerm(InnerQuantumState InnerState, spinflipTerm *T)
{
    QuantumState in = S.getQuantumState(QN,InnerState);	//real state	
    if ( S.n_i(in,T->bit[0]) || S.n_i(in,T->bit[1]) || !S.n_i(in,T->bit[2]) || !S.n_i(in,T->bit[3]))
	return;
    QuantumState diff1 = (1<<T->bit[0]) + (1<<T->bit[1]);
    QuantumState diff2 = (1<<T->bit[2]) + (1<<T->bit[3]);
    if ((diff2 > in + diff1) || ((diff2 <= in) && (S.getNumberOfStates() - diff1 <= in - diff2)))
	return;
    QuantumState out = in + diff1 - diff2;	
    if (out>in)
	return;
#warning The fact that hamiltonian is hermitian isn't taken into account. There may be slight perfomance issues.
    QuantumNumbers out_info = S.getStateInfo(out);
    if (out_info==(QuantumNumbers) QN) 
    {
	InnerQuantumState OutInnerState = S.getInnerState(out);
	H(OutInnerState,InnerState) += T->Value*measurefunc(in,out,T->bit[0],T->bit[1],T->bit[2],T->bit[3]);
    }
};

int HamiltonianPart::measurefunc(QuantumState state1, QuantumState state2, int i, int j, int k, int l)	//auxiliary function
{
    int flag=1, p=0;
    for (int m=0; m<IndexInfo.getIndexSize();m++)	//checking of "hopping" between state1 and state2
    {
	if ( (m==i) || (m==j) || (m==k) || (m==l) )
	{
	    if ( (S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
	    if ( (S.n_i(state2,j)==0) || (S.n_i(state1,j)==1) ) {flag=0;break;} else flag=1;
	    if ( (S.n_i(state2,k)==1) || (S.n_i(state1,k)==0) ) {flag=0;break;} else flag=1;
	    if ( (S.n_i(state2,l)==1) || (S.n_i(state1,l)==0) ) {flag=0;break;} else flag=1;
	}
	else
	{
	    if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
	}
    }
    if (flag==1)	//choose of true sign
    {
	for (int m=0; m<j; m++) p+=S.n_i(state2,m);
	for (int s=0; s<i; s++) p+=S.n_i(state2,s);
	for (int t=0; t<k; t++) p+=S.n_i(state1,t);
	for (int u=0; u<l; u++) p+=S.n_i(state1,u);
    }
    return (flag*(1-2*(p%2)));
}

int HamiltonianPart::checkhop(long int state1, long int state2, int i, int j)			//analog measurefunc
{
    int flag=1, p=0;
    for (int m=0; m<IndexInfo.getIndexSize(); m++)
    {
	if( (m==i) || (m==j) )
	{
	    if ((S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
	    if ((S.n_i(state2,j)==1) || (S.n_i(state1,j)==0) ) {flag=0;break;} else flag=1;
	}
	else
	{
	    if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
	}
    }
    if (flag==1)
    {
	for (int m=0;m<i;m++) p+=S.n_i(state2,m);
	for (int u=0;u<j;u++) p+=S.n_i(state1,u);
    }
    return (flag*(1-2*(p%2)));
}

void HamiltonianPart::add_hopping(RealMatrixType& HoppingMatrix)
{
    for (int i=0;i<HoppingMatrix.rows();i++)
	for (int j=0;j<HoppingMatrix.cols();j++)
	    if (i!=j) add_hopping(i,j,HoppingMatrix(i,j));
}

void HamiltonianPart::add_hopping(int i, int j, RealType t)
{
    for ( InnerQuantumState st1=0; st1<H.rows(); st1++)
    {
	QuantumState state1 = S.getQuantumState(QN,st1);	//real state1
	QuantumState difference = (i>j)?(1<<i)-(1<<j):(1<<j)-(1<<i);
	if (!((difference > S.getNumberOfStates()-state1 && i>j) || (i<j && difference > state1 )))
 	{
	    QuantumState state2=(i>j)?state1+difference:state1-difference;
	    if (S.getStateInfo(state1) == S.getStateInfo(state2))
	    {
    		int st2 = S.getInnerState(state2);
		if (st2>=0) H(st1,st2)+= t*checkhop(state1,state2,i,j);
	    }
	}
    }
}

void HamiltonianPart::prepare()
{
    size_t BlockSize = S.getQuantumStates(QN).size();

    H.resize(BlockSize,BlockSize);		//creation of H[i][j]=0 
    H.setZero();

    for(InnerQuantumState st=0; st<BlockSize; st++)
    {
	// loop over terms
	std::list<Term*>::const_iterator it1;
	for ( it1=IndexInfo.getTermsList().getTerms(2).begin() ; it1 != IndexInfo.getTermsList().getTerms(2).end(); ++it1 )
	{
	    if (( *it1)->type == "n") add_nTerm(st,(nTerm*) *it1);
	};
	
	std::list<Term*>::const_iterator it2;
	for ( it2=IndexInfo.getTermsList().getTerms(4).begin() ; it2 != IndexInfo.getTermsList().getTerms(4).end(); ++it2 )
	{
	    if ( (*it2)->type == "nn") add_nnTerm(st,(nnTerm*) *it2);
	    else 
	    if ( (*it2)->type == "spinflip") add_spinflipTerm(st,(spinflipTerm*) *it2);
	};
    }
	
    add_hopping(IndexInfo.getHoppingMatrix());
	
    H.triangularView<Eigen::Lower>() = H.triangularView<Eigen::Upper>().transpose();
}

RealType HamiltonianPart::getMatrixElement(InnerQuantumState m, InnerQuantumState n) const	//return  H(m,n)
{
    return H(m,n);
}

RealType HamiltonianPart::getEigenValue(InnerQuantumState state) const // return Eigenvalues(state)
{
    return Eigenvalues(state);
}

InnerQuantumState HamiltonianPart::getSize(void) const
{
    return H.rows();
}

QuantumNumbers HamiltonianPart::getQuantumNumbers(void) const
{
    return QN;
}

BlockNumber HamiltonianPart::getBlockNumber(void) const
{
    return S.getBlockNumber(QN);
}

// //other functions

void HamiltonianPart::diagonalize()		//method of diagonalization classificated part of Hamiltonian
{
    if (H.rows() == 1)
    {
	Eigenvalues = H;
	H(0,0) = 1;
    }
    if (H.rows() > 1)
    {
	Eigen::SelfAdjointEigenSolver<RealMatrixType> Solver(H,Eigen::ComputeEigenvectors);
	H = Solver.eigenvectors();
	Eigenvalues = Solver.eigenvalues();	// eigenvectors are ready
    }
}

void HamiltonianPart::print_to_screen()		// print part of Hamiltonian to screen
{
    std::cout << H << std::endl;
    std::cout << std::endl;
}

RealVectorType HamiltonianPart::getEigenState(InnerQuantumState state) const
{
    return H.col(state);
}

RealType HamiltonianPart::getMinimumEigenvalue() const
{
    return Eigenvalues.minCoeff();
}

bool HamiltonianPart::reduce(RealType ActualCutoff)
{
    InnerQuantumState counter=0;
    for (counter=0; (counter< (unsigned int)Eigenvalues.size() && Eigenvalues[counter]<=ActualCutoff); ++counter){};
    std::cout << "Left " << counter << " eigenvalues : " << std::endl;
    if (counter) 
	{std::cout << Eigenvalues.head(counter) << std::endl << "_________" << std::endl;
	Eigenvalues = Eigenvalues.head(counter);
	H = H.topLeftCorner(counter,counter);
	return true;
    }
    else return false;
}

void HamiltonianPart::save(H5::CommonFG* RootGroup) const
{
    HDF5Storage::saveRealVector(RootGroup,"V",Eigenvalues);
    HDF5Storage::saveRealMatrix(RootGroup,"H",H);

    // Save QN
    hsize_t Dim[1] = {QN.size()};
    H5::DataSpace DataSpace(1,Dim);
    H5::DataSet DataSet = RootGroup->createDataSet("QN",H5::PredType::NATIVE_SHORT,DataSpace);

    short* buf = new short[Dim[0]];
    std::copy(QN.begin(), QN.end(), buf);
    DataSet.write(buf,H5::PredType::NATIVE_SHORT);
    delete [] buf;
}
 
void HamiltonianPart::load(const H5::CommonFG* RootGroup)
{
    HDF5Storage::loadRealVector(RootGroup,"V",Eigenvalues);
    HDF5Storage::loadRealMatrix(RootGroup,"H",H);

    // Load QN
    H5::DataSet DataSet = RootGroup->openDataSet("QN");

    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HamiltonianPart::load()","Unexpected multidimentional dataspace."));

    short* buf = new short[DataSpace.getSimpleExtentNpoints()];
    DataSet.read(buf,H5::PredType::NATIVE_SHORT);

    if(! (QN == QuantumNumbers(buf[0],buf[1],buf[2]))) // FIXME!!!
	throw(H5::DataSetIException("HamiltonianPart::load()",
				    "Data in the storage is for another set of quantum numbers."));

    delete [] buf;
}

} // end of namespace Pomerol

