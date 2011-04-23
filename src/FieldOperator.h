#ifndef __INCLUDE_FIELDOPERATOR_H
#define __INCLUDE_FIELDOPERATOR_H

#include"Misc.h"
#include"ComputableObject.h"
#include"output.h"
#include"StatesClassification.h"
#include"Hamiltonian.h" 
#include"FieldOperatorPart.h"

typedef std::pair<BlockNumber,BlockNumber> BlockMapping;

class OperatorContainer : public ComputableObject
{
protected:
    StatesClassification &System;
    Hamiltonian &H;

    unsigned short bit;
    std::vector<FieldOperatorPart*> Data;
    std::map<unsigned int,BlockNumber> mapPartsFromRight;        // A map from non-zero parts to their BlockNumber indices
    std::map<unsigned int,BlockNumber> mapPartsFromLeft;        // A map from output index to input index, hence there is a unique transform
    std::map<unsigned int,BlockNumber> mapRightToLeftIndex;        // A map from output index to input index, hence there is a unique transform
    std::map<unsigned int,BlockNumber> mapLeftToRightIndex;        // A map from output index to input index, hence there is a unique transform
    std::list<BlockMapping> LeftRightIndices;
    unsigned int size;

    virtual BlockNumber mapsTo(BlockNumber RightIndex)=0;
    virtual    QuantumNumbers mapsTo(QuantumNumbers in)=0;

public:
    OperatorContainer(StatesClassification &System, Hamiltonian &H, int bit);

    FieldOperatorPart& getPartFromLeftIndex(BlockNumber in);
    FieldOperatorPart& getPartFromLeftIndex(QuantumNumbers in);
    FieldOperatorPart& getPartFromRightIndex(BlockNumber out);
    FieldOperatorPart& getPartFromRightIndex(QuantumNumbers out);
    BlockNumber getLeftIndex(BlockNumber RightIndex);
    BlockNumber getRightIndex(BlockNumber LeftIndex);
    std::list<BlockMapping>& getNonTrivialIndices();

    virtual void prepare()=0;
    
    void compute();
    void dump();
    void print_to_screen();
    unsigned short getIndex() const;
};

class CreationOperator;
class AnnihilationOperator;

class CreationOperator : public OperatorContainer
{
    friend class AnnihilationOperator;
    BlockNumber mapsTo(BlockNumber RightIndex);
    QuantumNumbers mapsTo(QuantumNumbers in);
public:
    AnnihilationOperator& transpose();
        void prepare();
    
    CreationOperator(StatesClassification &System, Hamiltonian &H, int bit);
};

class AnnihilationOperator : public OperatorContainer
{
    friend class CreationOperator;
    BlockNumber mapsTo(BlockNumber RightIndex);
    QuantumNumbers mapsTo(QuantumNumbers in);
public:
    CreationOperator& transpose();
        void prepare();
    
    AnnihilationOperator(StatesClassification &System, Hamiltonian &H, int bit);
};

#endif // endif :: #ifdef __INCLUDE_FIELDOPERATOR_H
