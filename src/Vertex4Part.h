#ifndef ____DEFINE_VERTEX4_PART____
#define ____DEFINE_VERTEX4_PART____

#include <list>

#include "config.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "FieldOperator.h"
#include "DensityMatrixPart.h"

class Vertex4Part {
       
public:
  
    struct Vertex4TermType1{
        ComplexType Residue;
        ComplexType Poles[3];
        
        Vertex4TermType1(RealType weight, RealType E1, RealType E2, RealType E3, RealType E4, Permutation3& Permutation);
        
        // Residue/((Frequency1 - Poles[0])*(Frequency2 - Poles[1])*(-Frequency3 - Poles[2]))
        ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const;
    };  

private:
  
    FieldOperatorPart& O1;
    FieldOperatorPart& O2;
    FieldOperatorPart& O3;
    CreationOperatorPart& CX4;
    
    HamiltonianPart& Hpart1;
    HamiltonianPart& Hpart2;
    HamiltonianPart& Hpart3;
    HamiltonianPart& Hpart4;
    
    DensityMatrixPart& DMpart1; 
    DensityMatrixPart& DMpart2;
    DensityMatrixPart& DMpart3;
    DensityMatrixPart& DMpart4;

    Permutation3 Permutation;
    
    std::list<Vertex4TermType1> TermsType1;
  
public:
    Vertex4Part(FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation);
 
    void compute(void);
    void computeFromRight(void);
    void compute13(void);
    void compute22(void);
    ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3);

};

#endif // endif :: #ifndef ____DEFINE_VERTEX4_PART____
