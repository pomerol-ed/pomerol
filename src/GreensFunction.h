#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____

#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "config.h"
#include "iniconfig.h"
#include "output.h"
#include "StatesClassification.h"
#include "FieldOperator.h"
#include "DensityMatrix.h"
#include "GreensFunctionPart.h"

class GreensFunction {
    
    std::list<GreensFunctionPart*> parts;
    StatesClassification& S;
    Hamiltonian& H;
    AnnihilationOperator& C;
    CreationOperator& CX;
    DensityMatrix& DM;
       
    output_handle green_path;
    
public:
    GreensFunction(StatesClassification& S, Hamiltonian& H,
                   AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM,
                   output_handle &OUT);
    ~GreensFunction();
    
    void prepare(void);
    void compute(void);
    
    unsigned short getBit(size_t Position) const;
    
    ComplexType operator()(long MatsubaraNum);
    
    //void dump(void);
    string getPath();
    void dumpMatsubara(unsigned short points);
    
    std::list<GreensFunctionPart::GreensTerm> getTerms(void) const;
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

