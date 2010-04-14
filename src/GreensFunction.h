#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____

#include <sstream>

#include "config.h"
#include "iniconfig.h"
#include "output.h"
#include "StatesClassification.h"
#include "FieldOperator.h"
#include "DensityMatrix.h"
#include "GreensFunctionPart.h"

class GreensFunction {
    
    StatesClassification& S;
    GreensFunctionPart** parts;
    
    output_handle green_path;
    
public:
    GreensFunction(StatesClassification& S, Hamiltonian& H,
                   AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM,
                   output_handle &OUT);
    ~GreensFunction();
    
    ComplexType operator()(ComplexType Frequency);
    
    //void dump(void);
    string getPath();
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

