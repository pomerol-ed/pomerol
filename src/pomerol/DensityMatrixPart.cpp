#include "pomerol/DensityMatrixPart.hpp"

namespace Pomerol {

DensityMatrixPart::DensityMatrixPart(const HamiltonianPart& H, RealType beta, RealType GroundEnergy) :
    Thermal(beta), H(H), GroundEnergy(GroundEnergy), weights(H.getSize())
{}

RealType DensityMatrixPart::computeUnnormalized()
{
    weights = exp(-beta*(H.getEigenValues().array() - GroundEnergy));
    return weights.sum();
}

void DensityMatrixPart::normalize(RealType Z)
{
    weights /= Z;
    Z_part /= Z;
}

RealType DensityMatrixPart::getAverageEnergy() const
{
    return weights.dot(H.getEigenValues());
}

RealType DensityMatrixPart::getWeight(InnerQuantumState s) const
{
    return weights(static_cast<Eigen::Index>(s));
}

void DensityMatrixPart::truncate(RealType Tolerance)
{
    Retained = false;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s) {
        if(weights(static_cast<Eigen::Index>(s)) > Tolerance) {
            Retained = true;
            break;
        }
    }
}

} // namespace Pomerol
