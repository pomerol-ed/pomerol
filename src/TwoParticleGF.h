#ifndef ____DEFINE_2PGF____
#define ____DEFINE_2PGF____

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
#include "TwoParticleGFPart.h"

Permutation3 getPermutation3(size_t p);
Permutation4 getPermutation4(size_t p);

class TwoParticleGF {

    std::list<TwoParticleGFPart*> parts;

    StatesClassification& S;
    Hamiltonian& H;
    AnnihilationOperator& C1;
    AnnihilationOperator& C2;
    CreationOperator& CX3;
    CreationOperator& CX4;
    DensityMatrix& DM;

    output_handle green_path;

    FieldOperatorPart& OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in);
    BlockNumber getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex);
    BlockNumber getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex); //!< return right index of an operator at current position for a current permutation

    TwoParticleGFPart::MatsubaraContainer *Storage;

public:
    TwoParticleGF(StatesClassification& S, Hamiltonian& H,
            AnnihilationOperator& C1, AnnihilationOperator& C2,
            CreationOperator& CX3, CreationOperator& CX4,
            DensityMatrix& DM,
            output_handle &OUT);
    ~TwoParticleGF();

    void prepare(void);
    void compute(long NumberOfMatsubaras);

    unsigned short getBit(size_t Position) const;
    RealType getBeta() const;

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    string getPath();
    //void dumpMatsubara(unsigned short points);

    static const Permutation3 permutations3[2];
    static const Permutation4 permutations4[24];
};

#endif // endif :: #ifndef ____DEFINE_2PGF____

