#ifndef ____DEFINE_HAMILTONIAN_PART____
#define ____DEFINE_HAMILTONIAN_PART____
#include "config.h"
#include "BitClassification.h"
#include "StatesClassification.h"

class HamiltonianPart {

	BitClassification &Formula;
	StatesClassification &S;

	RealType F_0, F_2, Us, mu, mus, t, ts,U,J;	//U = F_0 +4*F_2/25, J = 3*F_2/25
	int ** W1; int ** W2; int ** W3;		//matrix {Wn}

	long int N_state_m;				//number states of vector st
	
	QuantumNumbers hpart_id;

	RealMatrixType H;				//part of Hamilt
	RealVectorType V;				//vector of Eigen Values

	string ef_path;					// EigenFunctions output path handler
	string ev_path;					// EigenVectors output path handler


private:

	void putmatrix();
	void putHamilt();
	
	//multiorbital functions

	int hop;					//fit4a for inhopfuncW_3
	
    void add_diag(int st, RealType F_0, double F_2);					     //function creates diagonal elements of Hamilt
	void add_nondiag(int st1, int st2, RealType F_2);					     //function creates nondiagonal elements of Hamilt
	
	
	int measurefunc(long int state1, long int state2, int i, int j, int k, int l);	     // basic function for next two functions
    int inhopfuncW_2(long int state1, long int state2,int i, int j);	             // function checks probability adding nondiagonal elements
	int inhopfuncW_3(long int state1, long int state2, int i, int j, int *p);	     // similar previous function
	
	// s-orbital functions
	
	void add_U(int st, RealType Us);

	//chem. potentials
	
	void add_mu(int st, RealType mu);						//adds chem. potential on multiorbital
	void add_mus(int st, RealType mus);						//adds chem. potential on s-orbitals

	//functions describe hoppings
	
	int checkhop(long int state1, long int state2, int i, int j);		//check probability hopping between state1 and state2
	int hoppingfunc(long int state1, long int state2, int i);		//return result of hopping from "i" in state2 to state1
	void add_hopping_everywhere(int st1,int st2, RealType t, double ts);			//function adds to Hamilt hopping electron from "i"
	void add_hopping(RealMatrixType& HoppingMatrix);			//function adds to Hamilt hopping electron from "i"
	void add_hopping(int i,int j, RealType t);			//function adds to Hamilt hopping electron from "i"

public:

	HamiltonianPart(BitClassification &F_, StatesClassification &S_, QuantumNumbers id_) : Formula(F_),S(S_),hpart_id(id_){};
	
	void iniHamiltonianPart( RealType J_c, double U_c, double Us_c, double mu_c, double mus_c, 
	RealType t_c, double ts_c, const string &ev_path_, const string &ef_path_);				 //initialization HamiltonianPart
	
    int size(void);
	RealType reH(int m, int n);		//return H(m,n)
	RealType reV(int m);			//return V(m)

	void diagonalization();			    //method of process diagonalization
	QuantumNumbers id();				//return id of current hpart
	
	void dump();				//writing Eigen Values and Eigen Vectors in output file
	void print_to_screen();			//print to screen part of hamiltonian
	
};

						//structure of values rotated C or CX
#endif // endif :: #ifndef ____DEFINE_HAMILTONIAN_PART____
