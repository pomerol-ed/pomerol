#ifndef ____DEFINE_GREEN____
#define ____DEFINE_GREEN____

#include<list>

#include "config.h"
#include "iniconfig.h"
#include "output.h"
#include "FieldOperators.h"

class green {

    struct greenTerm{
        ComplexType Residue;
        ComplexType Pole;
        
        greenTerm(ComplexType Residue, ComplexType Pole);
        ComplexType operator()(ComplexType Frequency) const;
    };
    
    std::list<greenTerm> Terms;
    
    output_handle green_path;
    
    Hamiltonian &H;
    AnnihilationOperator& C;
    CreationOperator& CX;
    
    void precompute(RealType beta);
    
public:
    green(Hamiltonian& H, AnnihilationOperator& C, CreationOperator& CX, output_handle &OUT);
    ComplexType operator()(ComplexType Frequency);
    
    void dump(void);
    string path();
};

#endif // endif :: #ifndef ____DEFINE_GREEN____

