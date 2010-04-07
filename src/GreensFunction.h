#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____

#include "config.h"
#include "output.h"
#include "getStates.h"
#include "FieldOperators.h"

class GreensFunction {
    
    GreensFunctionPart** parts;
    
    output_handle green_path;
    
public:
    GreensFunction(Hamiltonian& H, AnnihilationOperator& C, CreationOperator& CX, output_handle &OUT);
    ~GreensFunction()
    
    ComplexType operator()(ComplexType Frequency);
    
    void dump(void);
    string path();
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

