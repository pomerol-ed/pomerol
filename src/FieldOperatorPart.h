#ifndef ____DEFINE_CCXPART____
#define ____DEFINE_CCXPART____
#include "config.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "output.h"
#include "HamiltonianPart.h"

struct valC {
	QuantumState n;					//number of line of matrix C or CX
	QuantumState m;					//number of column of matrix C or CX
	RealType C;				//value of rotated C or CX
	valC(QuantumState line, QuantumState column, RealType C_nm);			//initialization valC
};

						//class rotates matrixes C and CX 

template<class StorageType> class FieldOperatorPart {
protected:
	int i;
	StorageType elements;			//vector of notrivial elements of rotated matrix C

	StatesClassification &S;
	HamiltonianPart &h_from;
	HamiltonianPart &h_to;
	output_handle OUT;			//output path handler

    // basic functions
    virtual QuantumState retK(QuantumState L)=0;  
    virtual int mFunc(QuantumState state1, QuantumState state2, int i)=0;   //checks matrix element of an operator between state1 and state2
    virtual bool checkL(QuantumState L)=0;  //checks state L to be appropriate as a result of a creation/destruction operator

public:
  
	FieldOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from,	HamiltonianPart &h_to, output_handle &OUT);

	void compute();
	void dump();
	void print_to_screen();						//print to screen matrices UXCU UXCXU

	const string& path();						//output paths
    
    StorageType& value();
};

class AnnihilationOperatorPart : public FieldOperatorPart<RowMajorMatrixType>
{ 
    QuantumState retK(QuantumState L);	
    int mFunc(QuantumState state1, QuantumState state2, int i);
    bool checkL(QuantumState L);

public :
    AnnihilationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle &OUT);
};

class CreationOperatorPart : public FieldOperatorPart<ColMajorMatrixType>
{
    QuantumState retK(QuantumState L);	
    int mFunc(QuantumState state1, QuantumState state2, int i);	
    bool checkL(QuantumState L);
  
public :
    CreationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle &OUT);
};

//class FieldOperatorPart                               //rotates matrixes C and CX 
#include <sstream>
#include <fstream>
#include <iomanip>

template<class StorageType> FieldOperatorPart<StorageType>::FieldOperatorPart(
        int i, StatesClassification &S, HamiltonianPart &h_from,  HamiltonianPart &h_to, output_handle &OUT) : 
        i(i), S(S), h_from(h_from), h_to(h_to), OUT(OUT)
{};

template<class StorageType> StorageType& FieldOperatorPart<StorageType>::value()
{
    return elements;
}

// other functions

template<class StorageType> void FieldOperatorPart<StorageType>::print_to_screen()  //print to screen C and CX
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

template<class StorageType> void FieldOperatorPart<StorageType>::dump() //writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
    std::stringstream filename;
    filename << (*this).OUT.fullpath() << "//" << "C" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
    ofstream outCpart;
    outCpart.open(filename.str().c_str());
        
    outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements << endl;

    outCpart.close();
    cout << "The part of field operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

template<class StorageType> const string &FieldOperatorPart<StorageType>::path()
{
  static string str=(*this).OUT.fullpath(); return str;
}

template<class StorageType> void FieldOperatorPart<StorageType>::compute()
{
QuantumNumbers to   = h_to.id();
QuantumNumbers from = h_from.id();
elements.resize(S.clstates(to).size(),S.clstates(from).size());

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


#endif // endif :: #ifdef ____DEFINE_CCXPAIR____
