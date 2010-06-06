#ifndef ____DEFINE_DENSITY_MATRIX_PART____
#define ____DEFINE_DENSITY_MATRIX_PART____

#include "HamiltonianPart.h"

class DensityMatrixPart
{
    HamiltonianPart& hpart;
    RealType beta;
    RealType GroundEnergy;
  
    QuantumState partSize;
    RealVectorType weights;
    RealType Z_part;
  
public:
    DensityMatrixPart(HamiltonianPart& hpart, RealType beta, RealType GroundEnergy);
    void normalize(RealType Z);
    
    RealType compute(void);
    RealType weight(int m);
};

#endif // endif :: #ifndef ____DEFINE_DENSITY_MATRIX_PART____