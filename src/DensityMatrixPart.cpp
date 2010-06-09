#include "DensityMatrixPart.h"

DensityMatrixPart::DensityMatrixPart(HamiltonianPart& hpart, RealType beta, RealType GroundEnergy) :
    hpart(hpart), beta(beta), GroundEnergy(GroundEnergy)
{
    partSize = hpart.size();
    weights.resize(partSize);
}

RealType DensityMatrixPart::compute(void)
{
    Z_part = 0;
    for(QuantumState m = 0; m < partSize; ++m){
        weights(m) = exp(-beta*(hpart.reV(m)-GroundEnergy));
        Z_part += weights(m);
    }
  
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

RealType DensityMatrixPart::getBeta(void)
{
    return beta;
}