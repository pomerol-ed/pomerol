#ifndef ____DEFINE_GREENS_FUNCTION_PART____
#define ____DEFINE_GREENS_FUNCTION_PART____

#include <list>

#include "config.h"
#include "getStates.h"
#include "hpart.h"
#include "FieldOperators.h"
#include "DensityMatrixPart.h"

class GreensFunctionPart
{
   struct greenTerm{
        ComplexType Residue;
        ComplexType Pole;
        
        greenTerm(ComplexType Residue, ComplexType Pole);
        ComplexType operator()(ComplexType Frequency) const;
    };
    
    std::list<greenTerm> Terms;
  
public:
    GreensFunctionPart(AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                       getHpart& Hpart, DensityMatrixPart& DMpart);
 
    ComplexType operator()(RealType Frequency);
};

#endif // endif :: #ifndef ____DEFINE_GREENS_FUNCTION_PART____
