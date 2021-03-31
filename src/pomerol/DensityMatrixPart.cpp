#include "pomerol/DensityMatrixPart.h"

namespace Pomerol{

template<bool Complex>
DensityMatrixPart<Complex>::DensityMatrixPart(const StatesClassification<Complex> &S,
                                              const HamiltonianPart<Complex>& hpart,
                                              RealType beta,
                                              RealType GroundEnergy) :
    Thermal(beta), S(S), hpart(hpart), GroundEnergy(GroundEnergy), weights(hpart.getSize()), retained(true)
{}

template<bool Complex>
RealType DensityMatrixPart<Complex>::computeUnnormalized(void)
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

template<bool Complex>
void DensityMatrixPart<Complex>::normalize(RealType Z)
{
    weights /= Z;
    Z_part /= Z;
}

template<bool Complex>
RealType DensityMatrixPart<Complex>::getPartialZ(void) const
{
    return Z_part;
}

template<bool Complex>
RealType DensityMatrixPart<Complex>::getAverageEnergy(void) const
{
    RealType E=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        E += weights(s)*hpart.getEigenValue(s);
    }
    return E;
};

template<bool Complex>
RealType DensityMatrixPart<Complex>::getAverageOccupancy(void) const
{
    RealType n=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        VectorType<Complex> CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            n += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi).count()*
		    std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    };
    return n;
};

template<bool Complex>
RealType DensityMatrixPart<Complex>::getAverageOccupancy(ParticleIndex i) const
{
    RealType n=0.;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){
        VectorType<Complex> CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            n += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi).test(i)*
		    std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    };
    return n;
};

template<bool Complex>
RealType DensityMatrixPart<Complex>::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const
{
    RealType NN=0.;
    QuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s){ // s is an EigenState number
        VectorType<Complex> CurrentEigenState = hpart.getEigenState(s);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            NN += weights(s)*
		    S.getFockState(hpart.getBlockNumber(),fi)[i]*
		    S.getFockState(hpart.getBlockNumber(),fi)[j]*
		std::abs(CurrentEigenState(fi)*CurrentEigenState(fi));
    }
    return NN;
};

template<bool Complex>
RealType DensityMatrixPart<Complex>::getWeight(InnerQuantumState s) const
{
    return weights(s);
}

template<bool Complex>
void DensityMatrixPart<Complex>::truncate(RealType Tolerance)
{
    retained = false;
    InnerQuantumState partSize = weights.size();
    for(InnerQuantumState s = 0; s < partSize; ++s)
        if ( weights(s) > Tolerance ){
            retained = true;
            break;
        }
}

template<bool Complex>
bool DensityMatrixPart<Complex>::isRetained() const
{
    return retained;
}

template class DensityMatrixPart<false>;
template class DensityMatrixPart<true>;

} // end of namespace Pomerol
