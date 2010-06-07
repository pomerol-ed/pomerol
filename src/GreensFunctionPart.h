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
    HamiltonianPart& HpartInner;
    HamiltonianPart& HpartOuter;
    DensityMatrixPart& DMpartInner;
    DensityMatrixPart& DMpartOuter;
   
    AnnihilationOperatorPart& C;
    CreationOperatorPart& CX;
            
public:
    struct GreensTerm{
    ComplexType Residue;
        ComplexType Pole;
        
        GreensTerm(ComplexType Residue, ComplexType Pole);
        ComplexType operator()(ComplexType Frequency) const;
    };  
  
private:
    std::list<GreensTerm> Terms;
  
public:
    
    GreensFunctionPart(AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                       HamiltonianPart& HpartInner, HamiltonianPart& HpartOuter,
                       DensityMatrixPart& DMpartInner, DensityMatrixPart& DMpartOuter);
 
    void compute(void);
    ComplexType operator()(ComplexType Frequency) const;
    
    const std::list<GreensTerm>& getTerms(void) const;
};

std::ostream& operator<< (std::ostream& out, const GreensFunctionPart::GreensTerm& Term);

#endif // endif :: #ifndef ____DEFINE_GREENS_FUNCTION_PART____
