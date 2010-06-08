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
  
    enum ComputationMethod {ChasingIndices0, ChasingIndices1, ChasingIndices2};
    enum Var {Var1 = 0, Var2 = 1, Var3 = 2}; 
    
    struct Vertex4TermType1{
        ComplexType Residue;
        ComplexType Poles[3];
        
        Vertex4TermType1(RealType weight, RealType E1, RealType E2, RealType E3, RealType E4, Permutation3& Permutation);
        
        // Residue/((Frequency1 - Poles[0])*(Frequency2 - Poles[1])*(-Frequency3 - Poles[2]))
        ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const;
    };
    
    struct Vertex4TermType2{      
        Var z1, z2, z3;
        
        RealType Weight1PlusWeight2;
        RealType Weight1PlusWeight4;
        RealType BetaWeight2;
        
        RealType E4MinusE3;
        RealType E4MinusE2;
        RealType E2MinusE1;
        RealType E4MinusE1;
        
        bool IsDeltaE2E4Vanishes;
        
        Vertex4TermType2(RealType weight1, RealType weight2, RealType weight4, RealType beta, 
                         RealType E1, RealType E2, RealType E3, RealType E4, Permutation3& Permutation);
                         
        // 1) z2+z3 == 0 && E2 == E4
        //    1/((z3-E4MinusE3)*(z1-E2MinusE1))*( Weight1PlusWeight2/(z1-E2MinusE1) - BetaWeight2)
        // 2) otherwise:
        //    1/((z3-E4MinusE3)*(z2+z3-E4MinusE2))*(Weight1PlusWeight2/(z1-E2MinusE1) + Weight1PlusWeight4/(z1+z2+z3-E4MinusE1))
        ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const;
    };
    
    struct Vertex4TermType3{
        Var z1, z2, z3;
      
        RealType CoeffResonant;
        RealType CoeffNonResonant;
        
        RealType E4MinusE3;
        RealType E3MinusE2;
        RealType E3MinusE1;
        
        bool IsDeltaE2E3Vanishes;
        
        Vertex4TermType3(RealType weight1, RealType weight3, RealType beta, 
                         RealType E1, RealType E2, RealType E3, RealType E4, Permutation3& Permutation);
        
        // 1) z1+z2 == 0 && E2 == E3
        //    CoeffResonant/((z3-E4MinusE3)*(z2-E3MinusE2))
        // 2) otherwise:
        //    CoeffNonResonant/((z3-E4MinusE3)*(z2-E3MinusE2)*(z1+z2-E3MinusE1))
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
    
    void computeChasing0(void);
    void computeChasing1(void);
    void computeChasing2(void);
    
    void computeReduciblePart(void);
  
public:
    Vertex4Part(FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation);
 
    void compute(ComputationMethod method = ChasingIndices2);
    ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3);

};

#endif // endif :: #ifndef ____DEFINE_VERTEX4_PART____
