#ifndef ____DEFINE_GETSTATES____
#define ____DEFINE_GETSTATES____
#include "config.h"
#include "BitClassification.h"

struct BlockNumber {

	int number;

	BlockNumber(){};
	BlockNumber(int number_):number(number_){};
	operator int() const {return number;}
	BlockNumber& operator ++(int unused){number++; return *this;}

	bool isCorrect(){return number >= 0;}
	bool operator<(const BlockNumber& rhs) const ;
	bool operator==(const BlockNumber& rhs) const;
};

const BlockNumber ERROR_BLOCK_NUMBER = -1;

struct QuantumNumbers {			

	int Lz;					
	int N_up;
	int N_down;

	QuantumNumbers(int LZ, int N_UP, int N_DOWN);
	QuantumNumbers();
	friend std::ostream& operator<<(std::ostream& output, const QuantumNumbers& out);
	bool operator==(const QuantumNumbers &rhs)const ;
	bool operator<(const QuantumNumbers &rhs) const ;
};

const QuantumNumbers ERROR_QUANTUM_NUMBERS = QuantumNumbers(0,-1,-1);

typedef unsigned long int QuantumState;
typedef unsigned long int InnerQuantumState;

class StatesClassification {

	QuantumState N_state;			//2^N_bit number of states

	int N_bit;				//number bit of states
	int N_bit_m;				//2*(2*L(orbital) +1)

	vector<QuantumState> *** st;		//massive of vectors of states with Lz = "Lz", N_up = "N_up", N_down = "N_down" 
	int size;				//number of classificated vectors
	std::map<BlockNumber,QuantumNumbers> BlockToQuantum;
	std::map<QuantumNumbers,BlockNumber> QuantumToBlock;

	BlockNumber maximumBlockNumber_; 
	unsigned int BLOCKNUMBERLIMIT;
	BitClassification &Formula;

	void putstates();			//function gets all clasificated states with Lz = "Lz", N_up = "N_up", N_down = "N_down"
public:		
	StatesClassification(BitClassification& Formula):Formula(Formula) {};

	void iniStatesClassification();     				//inicialization StatesClassification

	const int N_b();							//return N_bit
	const int N_b_m();							//return N_bit_m
	const int L();								//return value of orbital moment
	const QuantumState N_st();							//return N_state

	const vector<QuantumState>& clstates( QuantumNumbers in );				//return st[in.Lz][in.N_up][in.N_down]
	const QuantumState cst( QuantumNumbers in, int m);					//return st[in.Lz][in.N_up][in.N_down][m]
	const InnerQuantumState getInnerState( QuantumState state);					//finds number of state @state in corresponding block 

	BlockNumber getBlockNumber(QuantumNumbers in);			//returns a number of Block which corresponds to given Quantum Numbers
	QuantumNumbers getBlockInfo(BlockNumber in);		//return Lz,N_up,N_down of number num
	BlockNumber NumberOfBlocks();						//return amount of non-trivial hamiltonian blocks
	QuantumNumbers getStateInfo(QuantumState in);		//return Lz,N_up,N_down of number num
	BlockNumber getBlockNumber(QuantumState in);			//returns a number of Block which corresponds to given Quantum Numbers

	inline int n_i(long int state, int i) { return ((state&(1<<i))>>i) ;};
	void getSiteInfo(int bit, int& lz, int& spin);

	bool checkQuantumNumbers(QuantumNumbers in);
};

#endif // endif :: #ifndef ____DEFINE_GETSTATES____
