#ifndef ____DEFINE_GETSTATES____
#define ____DEFINE_GETSTATES____
#include "config.h"

struct BlockNumber {

	int number;

	BlockNumber(){};
	BlockNumber(int number_):number(number_){};
	operator int(){return number;}
	BlockNumber& operator ++(int unused){number++; return *this;}

	bool isCorrect(){return number >= 0;}
};

const BlockNumber ERROR_BLOCK_NUMBER = -1;

struct QuantumNumbers {			

	int Lz;					
	int N_up;
	int N_down;

	QuantumNumbers(int LZ, int N_UP, int N_DOWN);
	QuantumNumbers();
	friend std::ostream& operator<<(std::ostream& output, const QuantumNumbers& out);
	bool operator==(QuantumNumbers &rhs){return (rhs.Lz == (*this).Lz && rhs.N_up == (*this).N_up && (*this).N_down == rhs.N_down );}
	bool operator==(QuantumNumbers rhs){return (rhs.Lz == (*this).Lz && rhs.N_up == (*this).N_up && (*this).N_down == rhs.N_down );}
};

const QuantumNumbers ERROR_QUANTUM_NUMBERS = QuantumNumbers(0,-1,-1);

typedef unsigned long int QuantumState;

class StatesClassification {

	QuantumState N_state;			//2^N_bit number of states

	int N_bit;				//number bit of states
	int N_bit_m;				//2*(2*L(orbital) +1)

	vector<QuantumState> *** st;		//massive of vectors of states with Lz = "Lz", N_up = "N_up", N_down = "N_down" 
	int size;				//number of classificated vectors
	vector<BlockNumber> num_bl;			//index of notrivial block,that appropriate defined Lz,N_up,N_down
	vector<QuantumNumbers> blockInfo;

	BlockNumber maximumBlockNumber_; 
	unsigned int BLOCKNUMBERLIMIT;

	void putstates();			//function gets all clasificated states with Lz = "Lz", N_up = "N_up", N_down = "N_down"
public:		
	StatesClassification() {}

	void iniStatesClassification(int N_BIT, int N_BIT_M);     				//inicialization StatesClassification

	const int N_b();							//return N_bit
	const int N_b_m();							//return N_bit_m
	const int L();								//return value of orbital moment
	const QuantumState N_st();							//return N_state

	const vector<QuantumState>& clstates( QuantumNumbers in );				//return st[in.Lz][in.N_up][in.N_down]
	const QuantumState cst( QuantumNumbers in, int m);					//return st[in.Lz][in.N_up][in.N_down][m]
	const long int inner_state( QuantumState state);					//finds number of state @state in corresponding block 

	BlockNumber getBlockNumber(QuantumNumbers in);			//returns a number of Block which corresponds to given Quantum Numbers
	QuantumNumbers getBlockInfo(BlockNumber in);		//return Lz,N_up,N_down of number num
	BlockNumber NumberOfBlocks();						//return amount of non-trivial hamiltonian blocks
	QuantumNumbers getStateInfo(int num);		//return Lz,N_up,N_down of number num

	int n_i(long int state,int i);	 	//checks if electron sit on "i" in "state"
	void getSiteInfo(int bit, int& lz, int& spin);
};

#endif // endif :: #ifndef ____DEFINE_GETSTATES____
