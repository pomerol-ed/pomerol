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

template<class PartType> 
OperatorContainer<PartType>::OperatorContainer(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit) : 
    System(System), H(H), OUT(OUT), bit(bit)
{
    size=0;
}

template<class PartType>
void OperatorContainer<PartType>::prepare()
{
  Data = new PartType* [System.NumberOfBlocks()];
  for (BlockNumber b=0;b<System.NumberOfBlocks();b++)
    {
      if (where(b).isCorrect()) 
      {
     Data[b]=new PartType(bit,System,H.part(b),H.part(where(b)),OUT);
         cout << "Entering " << operatorName << " Operator part " << System.getBlockInfo(b) << "->" << System.getBlockInfo(where(b)) << endl; 
         mapNontrivialParts[size]=b;
         mapLeftToRightPart[where(b)]=b;
         size++;
      }    
    }
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromRightIndex(BlockNumber in)
{
  return *Data[in];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromRightIndex(QuantumNumbers in)
{
  return *Data[System.getBlockNumber(in)];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromLeftIndex(BlockNumber in)
{
  return *Data[mapLeftToRightPart[in]];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromLeftIndex(QuantumNumbers in)
{
  return *Data[mapLeftToRightPart[System.getBlockNumber(in)]];
}

template<class PartType>
void OperatorContainer<PartType>::print_to_screen()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
        Data[mapNontrivialParts[b_in]]->print_to_screen();
  };
}

template<class PartType>
void OperatorContainer<PartType>::compute()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
        Data[mapNontrivialParts[b_in]]->compute();
  };
}

template<class PartType>
void OperatorContainer<PartType>::dump()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
        Data[mapNontrivialParts[b_in]]->dump();
  };
}

#endif
