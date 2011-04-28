/** \file src/DensityMatrixPart.cpp
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "DensityMatrixPart.h"

DensityMatrixPart::DensityMatrixPart(StatesClassification &S, HamiltonianPart& hpart, RealType beta, RealType GroundEnergy) :
    S(S), hpart(hpart), beta(beta), GroundEnergy(GroundEnergy)
{
    partSize = hpart.size();
    weights.resize(partSize);
}

RealType DensityMatrixPart::compute(void)
{
    Z_part = 0;
    for(QuantumState m = 0; m < partSize; ++m){
        // The non-normalized weight is <=1 for any state.
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

RealType DensityMatrixPart::getAverageEnergy()
{
    RealType E=0.;
    for(QuantumState m = 0; m < partSize; ++m){
        E += weights(m)*hpart.reV(m);
    }
    return E;
};

RealType DensityMatrixPart::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j)
{
    RealType NN=0.;
    for(InnerQuantumState m = 0; m < partSize; ++m){ // m is an EigenState number
        RealVectorType CurrentEigenState = hpart.getEigenState(m);
        for (InnerQuantumState fi=0; (long) fi < CurrentEigenState.size(); ++fi)
            NN += weights(m)*S.n_i(S.cst(hpart.id(),fi),i)*S.n_i(S.cst(hpart.id(),fi),j)*CurrentEigenState(fi)*CurrentEigenState(fi);
    }
    return NN;
};

/*
InnerQuantumState DensityMatrixPart::getMaximumTruncationState(RealType TruncationTolerance)
{
    InnerQuantumState m=0;
    DEBUG("");
    for (m=0; m<partSize && weights(m)>TruncationTolerance; ++m) DEBUG(hpart.reV(m) << " -> " << weights(m));
    DEBUG("m = " << m << " size = " << partSize << endl << endl);
    return m; // m>partSize-1?partSize-1:m;
}
*/

RealType DensityMatrixPart::weight(int m)
{
    return weights(m);
}

RealType DensityMatrixPart::getBeta(void)
{
    return beta;
}

#ifdef pomerolHDF5
void DensityMatrixPart::dumpIt(H5::CommonFG* FG) const
{
    Dumper::dumpReal(*FG,"beta",beta);
    Dumper::dumpReal(*FG,"GroundEnergy",GroundEnergy);
    Dumper::dumpReal(*FG,"Z_part",Z_part);
    Dumper::dumpRealVector(*FG,"weights",weights);
}

#endif // endif :: #ifdef pomerolHDF5
