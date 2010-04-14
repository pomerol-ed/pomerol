#include "DensityMatrixPart.h"

DensityMatrixPart::DensityMatrixPart(HamiltonianPart& hpart, RealType beta)
{
    QuantumState partSize = hpart.size();
  
    weights.resize(partSize);
    
    Z_part = 0;
    for(QuantumState m = 0; m < partSize; ++m){
        weights(m) = exp(-beta*hpart.reV(m));
        Z_part += weights(m);
    }
}

RealType DensityMatrixPart::getZ(void)
{
    return Z_part;
}

void DensityMatrixPart::normalize(RealType Z)
{
    weights /= Z;
    Z_part /= Z;
}

RealType DensityMatrixPart::weight(int m)
{
    return weights(m);
}