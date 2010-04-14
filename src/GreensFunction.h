#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____

#include "config.h"
#include "output.h"
#include "StatesClassification.h"
#include "FieldOperator.h"
#include "DensityMatrix.h"
#include "GreensFunctionPart.h"

class GreensFunction {
    
    GreensFunctionPart** parts;
    
    output_handle green_path;
    
public:
    GreensFunction(Hamiltonian& H, AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM, output_handle &OUT);
    ~GreensFunction();
    
    ComplexType operator()(ComplexType Frequency);
    
    void dump(void);
    string path();
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

