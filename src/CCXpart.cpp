#include "CCXpart.h"
#include <sstream>
#include <fstream>
#include <iomanip>

using std::stringstream;

//struct valC								//values rotated C or CX

valC::valC(QuantumState line, QuantumState column, RealType C_nm)	//inicialization valC
{	
	n=line;
	m=column;
	C=C_nm;
	
}

//class FieldOperatorPart								//rotates matrixes C and CX 

RealSparseMatrixType& FieldOperatorPart::value()						//return UXCU
{
	return elements;
}

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

// other functions

void FieldOperatorPart::print_to_screen()						//print to screen C and CX
{
	QuantumNumbers to   = h_to.id();
	QuantumNumbers from = h_from.id();
	for (int P=0; P<elements.outerSize(); ++P)
		for (RealSparseMatrixType::InnerIterator it(elements,P); it; ++it)
		{
				QuantumState N = S.clstates(to)[it.row()];
				QuantumState M = S.clstates(from)[it.col()];
				cout << N <<" " << M << " : " << it.value() << endl;
		};
}	

bool AnnihilationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==0);
}

bool CreationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==1);
}

void FieldOperatorPart::dump()							//writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
	stringstream filename;
	filename << (*this).OUT.fullpath() << "//" << "C" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
  	ofstream outCpart;
	outCpart.open(filename.str().c_str());
		
	outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements << endl;

	outCpart.close();
	cout << "The part of creation operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

const string &FieldOperatorPart::path()
{
  static string str=(*this).OUT.fullpath(); return str;
}

void FieldOperatorPart::compute()
{
QuantumNumbers to   = h_to.id();
QuantumNumbers from = h_from.id();
// elements = RealSparseMatrixType(S.clstates(to).size(),S.clstates(from).size());
elements.resize(S.clstates(to).size(),S.clstates(from).size());

for (std::vector<QuantumState>::const_iterator current_state=S.clstates(to).begin() ; current_state < S.clstates(to).end(); current_state++)
{

	QuantumState L=*current_state;
		
	if (checkL(L))
	{
		int K = retK(L);

		if( (mFunc(L,K,i)!= 0) )
		{			
								
			int l=S.inner_state(L), k=S.inner_state(K);				// l,k in part of Hamilt			
				
			for ( unsigned int n=0; n<S.clstates(to).size(); n++)
			{
				if(h_to.reH(l,n)!=0)
				{
					for (unsigned int m=0; m<S.clstates(from).size(); m++)
					{
						RealType C_nm = h_to.reH(l,n)*mFunc(L,K,i)*h_from.reH(k,m);
						if (fabs(C_nm)>MATRIX_ELEMENT_TOLERANCE)
						{	
							elements.coeffRef(n,m)+=C_nm;
						}
					}
				}
			}	
		}
	}	
};
  cout << elements << endl;
 elements.prune(MATRIX_ELEMENT_TOLERANCE);
  cout << elements << endl;
}
