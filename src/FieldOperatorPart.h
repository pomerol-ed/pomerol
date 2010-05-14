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

    RealType computeElement(const QuantumState row, const QuantumState col, const QuantumNumbers &from, const QuantumNumbers &to);
    
public:
  
	FieldOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from,	HamiltonianPart &h_to, output_handle OUT);

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
    AnnihilationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT);
};

class CreationOperatorPart : public FieldOperatorPart<ColMajorMatrixType>
{
    QuantumState retK(QuantumState L);	
    int mFunc(QuantumState state1, QuantumState state2, int i);	
    bool checkL(QuantumState L);
  
public :
    CreationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT);
};

#endif // endif :: #ifdef ____DEFINE_CCXPAIR____