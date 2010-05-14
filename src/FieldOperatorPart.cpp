#include "FieldOperatorPart.h"
#include "FieldOperatorPart.tmpl.h"

using std::stringstream;

//struct valC								//values rotated C or CX

valC::valC(QuantumState line, QuantumState column, RealType C_nm)	//inicialization valC
{
	n=line;
	m=column;
	C=C_nm;
}

// Functions of specialized classes

AnnihilationOperatorPart::AnnihilationOperatorPart(int i, StatesClassification &S, 
                                                   HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT) :
FieldOperatorPart<RowMajorMatrixType>(i,S,h_from,h_to,output_handle(OUT.path()+"/matrixC"))
{}

CreationOperatorPart::CreationOperatorPart(int i, StatesClassification &S, 
                                                  HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT) :
FieldOperatorPart<ColMajorMatrixType>(i,S,h_from,h_to,output_handle(OUT.path()+"/matrixCX"))
{}

template void FieldOperatorPart<RowMajorMatrixType>::dump();
template void FieldOperatorPart<ColMajorMatrixType>::dump();

template RowMajorMatrixType& FieldOperatorPart<RowMajorMatrixType>::value();
template ColMajorMatrixType& FieldOperatorPart<ColMajorMatrixType>::value();

QuantumState AnnihilationOperatorPart::retK(QuantumState L)							//return K for C
{	
	return( L + (1<<i) );
}

QuantumState CreationOperatorPart::retK(QuantumState L)							//return K for CX
{	

	return( L - (1<<i) );

}

int AnnihilationOperatorPart::mFunc(QuantumState state1, QuantumState state2, int i)
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b(); m++)
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

int CreationOperatorPart::mFunc(QuantumState state1, QuantumState state2, int i)
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b(); m++)
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
