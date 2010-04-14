#ifndef ____DEFINE_DENSITY_MATRIX_PART____
#define ____DEFINE_DENSITY_MATRIX_PART____

#include "HamiltonianPart.h"

class DensityMatrixPart
{
    RealVectorType weights;
    RealType Z_part;
  
public:
    DensityMatrixPart(HamiltonianPart& hpart, RealType beta);
    void normalize(RealType Z);
    
    RealType getZ(void);
    RealType weight(int m);
};

#endif // endif :: #ifndef ____DEFINE_DENSITY_MATRIX_PART____