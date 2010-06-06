#ifndef ____DEFINE_CCXPART____
#define ____DEFINE_CCXPART____
#include "config.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "output.h"
#include "HamiltonianPart.h"

class FieldOperatorPart {
protected:
	unsigned short i;
	SparseMatrixType elements;	//vector of notrivial elements of rotated matrix C

	StatesClassification &S;
	HamiltonianPart &h_from;
	HamiltonianPart &h_to;
	output_handle OUT;			//output path handler

    // basic functions
    virtual QuantumState retK(QuantumState L)=0;  
    virtual int mFunc(QuantumState state1, QuantumState state2, int i)=0;   //checks matrix element of an operator between state1 and state2
    virtual bool checkL(QuantumState L)=0;  //checks state L to be appropriate as a result of a creation/destruction operator
   
public:
  
	FieldOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from,	HamiltonianPart &h_to, output_handle OUT);

	void compute();
	void dump();
	void print_to_screen();						//print to screen matrices UXCU UXCXU

	const string& path();						//output paths
    
    SparseMatrixType& value();
};

class AnnihilationOperatorPart : public FieldOperatorPart
{ 
    QuantumState retK(QuantumState L);	
    int mFunc(QuantumState state1, QuantumState state2, int i);
    bool checkL(QuantumState L);

public :
    AnnihilationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT);
};

class CreationOperatorPart : public FieldOperatorPart
{
    QuantumState retK(QuantumState L);	
    int mFunc(QuantumState state1, QuantumState state2, int i);	
    bool checkL(QuantumState L);
  
public :
    CreationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to, output_handle OUT);
};

#endif // endif :: #ifdef ____DEFINE_CCXPAIR____