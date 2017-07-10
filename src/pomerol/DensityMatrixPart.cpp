#include "pomerol/DensityMatrixPart.h"

namespace Pomerol{
DensityMatrixPart::DensityMatrixPart(const StatesClassification &S, const HamiltonianPart& hpart, RealType beta, RealType GroundEnergy) :
    Thermal(beta), S(S), hpart(hpart), GroundEnergy(GroundEnergy), weights(hpart.getSize())
{}

RealType DensityMatrixPart::computeUnnormalized(void)
{
    Z_part = 0;
    QuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        // The non-normalized weight is <=1 for any state.
        weights(s) = exp(-beta*(hpart.getEigenValue(s)-GroundEnergy));
        Z_part += weights(s);
    }
    return Z_part;
}

void DensityMatrixPart::normalize(RealType Z)
{
    weights /= Z;
    Z_part /= Z;
}

RealType DensityMatrixPart::getPartialZ(void) const
{
    return Z_part;
}

RealType DensityMatrixPart::getAverageEnergy(void) const
{
    RealType E=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        E += weights(s)*hpart.getEigenValue(s);
    }
    return E;
};

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

RealType DensityMatrixPart::getWeight(InnerQuantumState s) const
{
    return weights(s);
}

} // end of namespace Pomerol
