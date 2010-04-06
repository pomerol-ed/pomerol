#ifndef ____DEFINE_CCXPART____
#define ____DEFINE_CCXPART____
#include "config.h"
#include "getStates.h"
#include "Hamiltonian.h"
#include "output.h"
#include "hpart.h"

struct valC {
	QuantumState n;					//number of line of matrix C or CX
	QuantumState m;					//number of column of matrix C or CX
	RealType C;				//value of rotated C or CX
	valC(QuantumState line, QuantumState column, RealType C_nm);			//initialization valC
};

						//class rotates matrixes C and CX 

class FieldOperatorPart {
protected:
	int i;
	RealSparseMatrixType elements;			//vector of notrivial elements of rotated matrix C

	getStates &S;
	getHpart &h_from;
	getHpart &h_to;
	output_handle OUT;			//output path handler

	// basic functions
	virtual QuantumState retK(QuantumState L)=0;						
	virtual int mFunc(QuantumState state1, QuantumState state2, int i)=0;		//checks matrix element of an operator between state1 and state2	
	virtual bool checkL(QuantumState L)=0;						//checks state L to be appropriate as a result of a creation/destruction operator

public:
	FieldOperatorPart(
		int i_, 
		getStates &S_, 
		getHpart &h_from_,
		getHpart &h_to_, 
		output_handle &OUT_
	       ):i(i_), S(S_), h_from(h_from_),h_to(h_to_), OUT(OUT_) {};

  	void compute();
	void dump();
	void print_to_screen();						//print to screen matrices UXCU UXCXU

	const string& path();							//output paths

  	RealSparseMatrixType &value();
};

class AnnihilationOperatorPart : public FieldOperatorPart
{
  QuantumState retK(QuantumState L);	
  int mFunc(QuantumState state1, QuantumState state2, int i);
  bool checkL(QuantumState L);
  public :
  AnnihilationOperatorPart(
  		  int i_,
  		  getStates &S_, 
                  getHpart &h_from_, 
		  getHpart &h_to_, 
		  output_handle &OUT_
		):FieldOperatorPart(i_,S_,h_from_,h_to_,OUT_){OUT=output_handle(OUT_.path()+"//matrixC");};
	    //   ):i(i_), S(S_), h_from(h_from_),h_to(h_to_) {OUT=output_handle(OUT_.path()+"//matrixC");}
};

class CreationOperatorPart : public FieldOperatorPart
{
  QuantumState retK(QuantumState L);	
  int mFunc(QuantumState state1, QuantumState state2, int i);	
  bool checkL(QuantumState L);
  public :
  CreationOperatorPart(
  		   int i_,
  		   getStates &S_, 
		   getHpart &h_from_, 
		   getHpart &h_to_,
		   output_handle &OUT_
		 ):FieldOperatorPart(i_,S_,h_from_,h_to_,OUT_){OUT=output_handle(OUT_.path()+"//matrixCX");};
};

#endif // endif :: #ifdef ____DEFINE_CCXPAIR____
