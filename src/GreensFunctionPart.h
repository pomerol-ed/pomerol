#ifndef ____DEFINE_GREENS_FUNCTION_PART____
#define ____DEFINE_GREENS_FUNCTION_PART____

#include <list>

#include "config.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "FieldOperator.h"
#include "DensityMatrixPart.h"

class GreensFunctionPart
{
   HamiltonianPart& Hpart;
   AnnihilationOperatorPart& C;
   CreationOperatorPart& CX;
   DensityMatrixPart& DMpart;
  
   struct greenTerm{
        ComplexType Residue;
        ComplexType Pole;
        
        greenTerm(ComplexType Residue, ComplexType Pole);
        ComplexType operator()(ComplexType Frequency) const;
    };
    
    std::list<greenTerm> Terms;
  
public:
    GreensFunctionPart(AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                       HamiltonianPart& Hpart, DensityMatrixPart& DMpart);
 
    void compute(void);
    ComplexType operator()(ComplexType Frequency);
};

#endif // endif :: #ifndef ____DEFINE_GREENS_FUNCTION_PART____
