#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____
#include "config.h"
#include "getStates.h"
#include "hpart.h"
#include "CCXpair.h"
#include "output.h"

class green {

	int i;					//number of M_sigma of C
	int j;					//number of M_sigma of CX
		
	ComplexType * G_ij;			//Green function

	output_handle green_path;		//Green function output path handler


	getStates &S;
	Hamiltonian &H;
	matrixs &operatorsCiCXj;
	output_handle &OUT;
	
private:

	void building();			//basic function - builds GF from matrixS and input

public:

	green(int I, int J, getStates &S_, Hamiltonian &H_, matrixs &operatorsCiCXj_, output_handle &OUT_);			//initialization class

	//other functions
	
	void dump();				//print to file Greens Function
	string path();				//returns path to Greens Function
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

