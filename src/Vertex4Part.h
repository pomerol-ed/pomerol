#ifndef ____DEFINE_VERTEX4_PART____
#define ____DEFINE_VERTEX4_PART____

#include <list>

#include "config.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "FieldOperator.h"
#include "DensityMatrixPart.h"

class Vertex4Part {
    
    AnnihilationOperatorPart& C0;
    AnnihilationOperatorPart& C1;
    CreationOperatorPart& CX2;
    CreationOperatorPart& CX3;
    
    HamiltonianPart& HpartOuter;
    HamiltonianPart& Hpart23;
    HamiltonianPart& Hpart12;
    HamiltonianPart& Hpart01;
    
    DensityMatrixPart& DMpartOuter; 
    DensityMatrixPart& DMpart23;
    DensityMatrixPart& DMpart12;
    DensityMatrixPart& DMpart01;
    
    size_t PermutationNumber;
    
    // TODO
    
public:
    Vertex4Part(AnnihilationOperatorPart& C0, AnnihilationOperatorPart& C1,
                CreationOperatorPart& CX2, CreationOperatorPart& CX3,
                HamiltonianPart& HpartOuter, HamiltonianPart& Hpart23,
                HamiltonianPart& Hpart12, HamiltonianPart& Hpart01,
                DensityMatrixPart& DMpartOuter, DensityMatrixPart& DMpart23,
                DensityMatrixPart& DMpart12, DensityMatrixPart& DMpart01,
                size_t PermutationNumber);
 
    void compute(void);
    ComplexType operator()(ComplexType Frequency0, ComplexType Frequency1, ComplexType Frequency2);

};

#endif // endif :: #ifndef ____DEFINE_VERTEX4_PART____