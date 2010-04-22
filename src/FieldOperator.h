#ifndef __OPERATOR_CONTAINERS__
#define __OPERATOR_CONTAINERS__

#include "config.h"
#include "output.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h" 
#include "FieldOperatorPart.h"

class OperatorContainer
{
protected:
	StatesClassification &System;
	Hamiltonian &H;
	output_handle OUT;

	int bit;
	FieldOperatorPart **Data;
	std::map<unsigned int,BlockNumber> mapNontrivialParts;		// A map from non-zero parts to their BlockNumber indices
	std::map<unsigned int,BlockNumber> mapLeftToRightPart;		// A map from output index to input index, hence there is a unique transform
	unsigned int size;

	virtual	BlockNumber where(BlockNumber in)=0;
	virtual	QuantumNumbers where(QuantumNumbers in)=0;

public:
	OperatorContainer(StatesClassification &System_, Hamiltonian &H_, output_handle &OUT_, int bit_):System(System_),H(H_),OUT(OUT_),bit(bit_){size=0;};	


	FieldOperatorPart& getPartLeftIndex(BlockNumber in);
	FieldOperatorPart& getPartLeftIndex(QuantumNumbers in);
	FieldOperatorPart& getPartRightIndex(BlockNumber out);
	FieldOperatorPart& getPartRightIndex(QuantumNumbers out);

	void compute();
	void dump();
	void print_to_screen();
};

class CreationOperator : public OperatorContainer
{
public:
	void prepare();
	BlockNumber where(BlockNumber in);
	QuantumNumbers where(QuantumNumbers in);

	CreationOperator(StatesClassification &System_, Hamiltonian &H_, output_handle &OUT_, int bit_):OperatorContainer(System_,H_,OUT_,bit_){};
};

class AnnihilationOperator : public OperatorContainer
{
public:
	void prepare();
	BlockNumber where(BlockNumber in);
	QuantumNumbers where(QuantumNumbers in);
	AnnihilationOperator(StatesClassification &System_, Hamiltonian &H_, output_handle &OUT_, int bit_):OperatorContainer(System_,H_,OUT_,bit_){};
};

#endif
