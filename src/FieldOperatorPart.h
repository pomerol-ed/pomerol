#ifndef __INCLUDE_FIELDOPERATORPART_H
#define __INCLUDE_FIELDOPERATORPART_H
#include"Misc.h"
#include"StatesClassification.h"
#include"Hamiltonian.h"
#include"output.h"

class FieldOperatorPart {
protected:
    unsigned short i;
    RowMajorMatrixType elementsRowMajor;    
    ColMajorMatrixType elementsColMajor;    

    StatesClassification &S;
    HamiltonianPart &h_from;
    HamiltonianPart &h_to;

    // basic functions
    virtual QuantumState retK(QuantumState L)=0;  
    virtual int mFunc(QuantumState state1, QuantumState state2, int i)=0;   //checks matrix element of an operator between state1 and state2
    virtual bool checkL(QuantumState L)=0;  //checks state L to be appropriate as a result of a creation/destruction operator
    
    static const RealType MatrixElementTolerance = 1e-8;
   
public:
  
    FieldOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from,    HamiltonianPart &h_to);

    void compute();
    void dump();
    void print_to_screen();                        //print to screen matrices UXCU UXCXU

    RowMajorMatrixType& getRowMajorValue();
    ColMajorMatrixType& getColMajorValue();
    BlockNumber getLeftIndex();
    BlockNumber getRightIndex();
};

class AnnihilationOperatorPart;
class CreationOperatorPart;

class AnnihilationOperatorPart : public FieldOperatorPart
{ 
    QuantumState retK(QuantumState L);    
    int mFunc(QuantumState state1, QuantumState state2, int i);
    bool checkL(QuantumState L);
    friend class CreationOperatorPart;

public :
    AnnihilationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to);
    CreationOperatorPart& transpose();
};

class CreationOperatorPart : public FieldOperatorPart
{
    QuantumState retK(QuantumState L);    
    int mFunc(QuantumState state1, QuantumState state2, int i);    
    bool checkL(QuantumState L);
    friend class AnnihilationOperatorPart;
  
public :
    CreationOperatorPart(int i, StatesClassification &S, HamiltonianPart &h_from, HamiltonianPart &h_to);
    AnnihilationOperatorPart& transpose();
};

#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H
