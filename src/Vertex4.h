#ifndef ____DEFINE_VERTEX4____
#define ____DEFINE_VERTEX4____

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
#include "Vertex4Part.h"

Permutation3 getPermutation3(size_t p);
Permutation4 getPermutation4(size_t p);

class Vertex4 {
    
    std::list<Vertex4Part*> parts;
    
    StatesClassification& S;
    Hamiltonian& H;
    AnnihilationOperator& C1;
    AnnihilationOperator& C2;
    CreationOperator& CX3;
    CreationOperator& CX4;
    DensityMatrix& DM;
    
    output_handle green_path;
    
    FieldOperatorPart& OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in);
    BlockNumber OperatorAtPositionMapsTo(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in);
    
public:
    Vertex4(StatesClassification& S, Hamiltonian& H,
            AnnihilationOperator& C1, AnnihilationOperator& C2,
            CreationOperator& CX3, CreationOperator& CX4,
            DensityMatrix& DM,
            output_handle &OUT);
    ~Vertex4();
    
    void prepare(void);
    void compute(void);
    
    ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3);
    
    string getPath();
    //void dumpMatsubara(unsigned short points);
    
    static const Permutation3 permutations3[2];
    static const Permutation4 permutations4[24];
};

#endif // endif :: #ifndef ____DEFINE_VERTEX4____

