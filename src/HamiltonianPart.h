#ifndef ____DEFINE_HAMILTONIAN_PART____
#define ____DEFINE_HAMILTONIAN_PART____
#include "config.h"
#include "BitClassification.h"
#include "StatesClassification.h"

class HamiltonianPart {

	BitClassification &Formula;
	StatesClassification &S;

	RealType mu, mus;	

	QuantumState N_state_m;				//number states of vector st
	
	QuantumNumbers hpart_id;

	RealMatrixType H;				//part of Hamilt
	RealVectorType V;				//vector of Eigen Values

	string ev_path;					// EigenVectors output path handler
	string ef_path;					// EigenFunctions output path handler


private:

	void add_nTerm(InnerQuantumState st, nTerm *N);

	void add_nnTerm(InnerQuantumState st, nnTerm *T);
	void add_spinflipTerm(InnerQuantumState st, spinflipTerm *T);
	
	int measurefunc(QuantumState state1, QuantumState state2, int i, int j, int k, int l);	     // basic function for next two functions
	
	// s-orbital functions

	//chem. potentials
	
	void add_mu(int st, RealType mu);						//adds chem. potential on multiorbital
	void add_mus(int st, RealType mus);						//adds chem. potential on s-orbitals

	//functions describe hoppings
	
	int checkhop(long int state1, long int state2, int i, int j);		//check probability hopping between state1 and state2
	void add_hopping(RealMatrixType& HoppingMatrix);			//function adds to Hamilt hopping electron from "i"
	void add_hopping(int i,int j, RealType t);			//function adds to Hamilt hopping electron from "i"

public:

	HamiltonianPart(BitClassification &F, StatesClassification &S, QuantumNumbers id, const string &ev_path, const string &ef_path) : Formula(F),S(S),hpart_id(id),ev_path(ev_path),ef_path(ef_path){};
	
	void enter( double mu_c, double mus_c);
	
    	InnerQuantumState size(void);
	RealType reH(int m, int n);		//return H(m,n)
	RealType reV(int m);			//return V(m)

	void diagonalization();			    //method of process diagonalization
	QuantumNumbers id();				//return id of current hpart
	
	void dump();				//writing Eigen Values and Eigen Vectors in output file
	void print_to_screen();			//print to screen part of hamiltonian
	
};

						//structure of values rotated C or CX
#endif // endif :: #ifndef ____DEFINE_HAMILTONIAN_PART____
