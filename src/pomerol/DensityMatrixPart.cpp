#include "pomerol/DensityMatrixPart.h"

namespace Pomerol{

DensityMatrixPart::DensityMatrixPart(const StatesClassification &S, const Hamiltonian& H, BlockNumber Block, RealType beta, RealType GroundEnergy) :
    Thermal(beta), S(S), H(H), Block(Block), GroundEnergy(GroundEnergy), weights(H.getBlockSize(Block)), retained(true)
{}

RealType DensityMatrixPart::computeUnnormalized(void)
{
    weights = exp(-beta*(H.getEigenValues(Block).array() - GroundEnergy));
    return weights.sum();
}

void DensityMatrixPart::normalize(RealType Z)
{
    weights /= Z;
    Z_part /= Z;
}

RealType DensityMatrixPart::getAverageEnergy() const
{
    return weights.dot(H.getEigenValues(Block));
}

// TODO
/*
RealType DensityMatrixPart::getAverageOccupancy(void) const
{
    RealType n=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        VectorType CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            n += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi).count()*
		    std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    };
    return n;
};

RealType DensityMatrixPart::getAverageOccupancy(ParticleIndex i) const
{
    RealType n=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        VectorType CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            n += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi).test(i)*
		    std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    };
    return n;
};

RealType DensityMatrixPart::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const
{
    RealType NN=0.;
    QuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){ // s is an EigenState number
        VectorType CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            NN += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi)[i]*
		    S.getFockState(hpart.getBlockNumber(),fi)[j]*
		std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    }
    return NN;
};
*/

RealType DensityMatrixPart::getWeight(InnerQuantumState s) const
{
    return weights(s);
}

void DensityMatrixPart::truncate(RealType Tolerance)
{
    retained = false;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s)
        if ( weights(s) > Tolerance ){
            retained = true;
            break;
        }
}

bool DensityMatrixPart::isRetained() const
{
    return retained;
}

} // end of namespace Pomerol
