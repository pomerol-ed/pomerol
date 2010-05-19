#ifndef __OPERATOR_CONTAINERS__
#define __OPERATOR_CONTAINERS__

#include "config.h"
#include "output.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h" 
#include "FieldOperatorPart.h"

template<class PartType> class OperatorContainer
{
protected:
    std::string operatorName;
  
	StatesClassification &System;
	Hamiltonian &H;
	output_handle OUT;

	int bit;
	PartType **Data;
	std::map<unsigned int,BlockNumber> mapNontrivialParts;		// A map from non-zero parts to their BlockNumber indices
	std::map<unsigned int,BlockNumber> mapLeftToRightPart;		// A map from output index to input index, hence there is a unique transform
	unsigned int size;

	virtual	BlockNumber where(BlockNumber in)=0;
	virtual	QuantumNumbers where(QuantumNumbers in)=0;

public:
	OperatorContainer(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit);

	PartType& getPartFromLeftIndex(BlockNumber in);
	PartType& getPartFromLeftIndex(QuantumNumbers in);
	PartType& getPartFromRightIndex(BlockNumber out);
	PartType& getPartFromRightIndex(QuantumNumbers out);

    void prepare();
    
	void compute();
	void dump();
	void print_to_screen();
};

class CreationOperator : public OperatorContainer<CreationOperatorPart>
{
public:
	BlockNumber where(BlockNumber in);
	QuantumNumbers where(QuantumNumbers in);

	CreationOperator(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit);
};

class AnnihilationOperator : public OperatorContainer<AnnihilationOperatorPart>
{
public:
	BlockNumber where(BlockNumber in);
	QuantumNumbers where(QuantumNumbers in);
	AnnihilationOperator(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit);
};

#include "FieldOperator.tmpl.h"

#endif
