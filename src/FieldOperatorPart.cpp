#include "FieldOperatorPart.h"
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

// class FieldOperatorPart
const string &FieldOperatorPart::path()
{
  static string str=OUT.fullpath(); return str;
}

// class AnnihilationOperatorPart

QuantumState AnnihilationOperatorPart::retK(QuantumState L)                         //return K for C
{   
    return( L + (1<<i) );
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

bool AnnihilationOperatorPart::checkL(QuantumState L)
{
    return (S.n_i(L,i)==0);
}

RowMajorMatrixType& AnnihilationOperatorPart::value()
{
    return elements;
}

void AnnihilationOperatorPart::print_to_screen()                        //print to screen C
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    for (int P=0; P<elements.outerSize(); ++P)
        for (RowMajorMatrixType::InnerIterator it(elements,P); it; ++it)
        {
                QuantumState N = S.clstates(to)[it.row()];
                QuantumState M = S.clstates(from)[it.col()];
                cout << N << " " << M << " : " << it.value() << endl;
        };
}

void AnnihilationOperatorPart::dump()           //writing AnnihilationOperatorPart C[M_sigma] in output file
{
    stringstream filename;
    filename << OUT.fullpath() << "/" << "C" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";

    ofstream outCpart(filename.str().c_str());
    outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements << endl;
    outCpart.close();

    cout << "The part of annihilation operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

void AnnihilationOperatorPart::compute()
{
QuantumNumbers to   = h_to.id();
QuantumNumbers from = h_from.id();
// elelements.resize(S.clstates(to).size(),S.clstates(from).size());
elements = RealSparseMatrixType(S.clstates(to).size(),S.clstates(from).size());

for (std::vector<QuantumState>::const_iterator current_state=S.clstates(to).begin() ; current_state < S.clstates(to).end(); current_state++)
{

    QuantumState L=*current_state;
        
    if (checkL(L))
    {
        int K = retK(L);

        if( (mFunc(L,K,i)!= 0) )
        {           
                                
            int l=S.inner_state(L), k=S.inner_state(K);             // l,k in part of Hamilt            
                
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
    elements.prune(MATRIX_ELEMENT_TOLERANCE);
}

// class CreationOperatorPart

QuantumState CreationOperatorPart::retK(QuantumState L)                         //return K for CX
{   
    return( L - (1<<i) );
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

bool CreationOperatorPart::checkL(QuantumState L)
{
    return (S.n_i(L,i)==1);
}

ColMajorMatrixType& CreationOperatorPart::value()                        //return UXCXU
{
    return elements;
}

void CreationOperatorPart::print_to_screen()						//print to screen CX
{
	QuantumNumbers to   = h_to.id();
	QuantumNumbers from = h_from.id();
	for (int P=0; P<elements.outerSize(); ++P)
		for (ColMajorMatrixType::InnerIterator it(elements,P); it; ++it)
		{
				QuantumState N = S.clstates(to)[it.row()];
				QuantumState M = S.clstates(from)[it.col()];
				cout << N <<" " << M << " : " << it.value() << endl;
		};
}	

void CreationOperatorPart::dump()			//writing CreationOperatorPart CX[M_sigma] t output file
{
	stringstream filename;
	filename << OUT.fullpath() << "/" << "CX" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";

    ofstream outCXpart(filename.str().c_str());
	outCXpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements << endl;
	outCXpart.close();

    cout << "The part of creation operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

void CreationOperatorPart::compute()
{
QuantumNumbers to   = h_to.id();
QuantumNumbers from = h_from.id();
// elelements.resize(S.clstates(to).size(),S.clstates(from).size());
elements = RealSparseMatrixType(S.clstates(to).size(),S.clstates(from).size());

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
    elements.prune(MATRIX_ELEMENT_TOLERANCE);
}