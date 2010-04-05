#ifndef ____DEFINE_CCXPAIR____
#define ____DEFINE_CCXPAIR____
#include "config.h"
#include "getStates.h"
#include "Hamiltonian.h"
#include "output.h"
#include "hpart.h"

struct valC {

	QuantumState 	n;					//number of line of matrix C or CX
	QuantumState    m;					//number of column of matrix C or CX
	
	RealType C;				//value of rotated C or CX
	
	valC& operator+=(const valC &rhs); 
	valC(int line, int column, RealType C_nm);			//inicialization valC
	
};

						//class rotates matrixes C and CX 

class matrixs {

	 					//M_sigma - "place" in state, where act C and CX
	int i;					//number of M_sigma of C
	int j;					//number of M_sigma of CX

	vector<valC> UXCU;			//vector of notrivial elements of rotated matrix C
	vector<valC> UXCXU;			//vector of notrivial elements of rotated matrix CX

	vector<valC> ** uxcu;			//classificated vector UXCU
	vector<valC> ** uxcxu;			//classificated vector UXCXU

	output_handle matrixC_path;		// matrixC output path handler
	output_handle matrixCX_path;		// matrixCX output path handler
	getStates &S;
	Hamiltonian &H;

	output_handle &OUT;


private:

	// basic functions
	
//	const vector<long int>& reVect( int st );			//return st[Lz][N_up][N_down]
//	getHpart& reHpart( int st );					//return getHpart[Lz][N_up][N_down]

	int retKforC(int L);						//return number of column (K) for C
	int retKforCX(int L);						//return number of column (K) for CX

	int mFuncC(long int state1, long int state2, int i);		//analog other measurefunc functions for operator C
	int mFuncCX(long int state1, long int state2, int i);		//analog other measurefunc functions for operator CX

public:
	matrixs(getStates &S_, Hamiltonian &H_, output_handle &OUT_):S(S_), H(H_), OUT(OUT_) {};
	
	void inimatrixs(int I, int J);					//inicialization class
	
	void putMatrXC();						//creation UXCU 
	void putMatrXCX();						//creation UXCXU
	
	vector<valC>& reVecC();						//return UXCU(m,n)
	vector<valC>& reVecCX();					//return UXCXU(m,n)

	//other functions
	
	void print_to_screen();						//print to screen matrices UXCU UXCXU
	void dump();							//print to file matrices UXCU UXCXU
	string path();							//output paths
};


#endif // endif :: #ifdef ____DEFINE_CCXPAIR____
