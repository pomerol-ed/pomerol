#ifndef ____DEFINE_DENSITY_MATRIX____
#define ____DEFINE_DENSITY_MATRIX____

#include <vector>

#include "StatesClassification.h"
#include "BitClassification.h"
#include "Hamiltonian.h"
#include "DensityMatrixPart.h"

class DensityMatrix
{
    StatesClassification& S;
    std::vector<DensityMatrixPart> parts;

public:
    DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta);
    
    RealType operator()( QuantumState &state );
};

#endif // endif :: #ifndef ____DEFINE_DENSITY_MATRIX____