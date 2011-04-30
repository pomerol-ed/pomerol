/** \file src/DensityMatrixPart.cpp
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "DensityMatrixPart.h"

DensityMatrixPart::DensityMatrixPart(StatesClassification &S, HamiltonianPart& hpart, RealType beta, RealType GroundEnergy) :
    Thermal(beta), S(S), hpart(hpart), GroundEnergy(GroundEnergy), weights(hpart.size())
{}

RealType DensityMatrixPart::compute(void)
{
    Z_part = 0;
    QuantumState partSize = weights.size();
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
    QuantumState partSize = weights.size();
    for(QuantumState m = 0; m < partSize; ++m){
        E += weights(m)*hpart.reV(m);
    }
    return E;
};

RealType DensityMatrixPart::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j)
{
    RealType NN=0.;
    QuantumState partSize = weights.size();
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

void DensityMatrixPart::save(H5::CommonFG* RootGroup, HDF5Storage const* const Storage) const
{
    Storage->saveReal(RootGroup,"beta",beta);
    Storage->saveReal(RootGroup,"GroundEnergy",GroundEnergy);
    Storage->saveReal(RootGroup,"Z_part",Z_part);
    Storage->saveRealVector(RootGroup,"weights",weights);
}

void DensityMatrixPart::load(const H5::CommonFG* RootGroup, HDF5Storage const* const Storage)
{
    RealType newBeta = Storage->loadReal(RootGroup,"beta");
    if(newBeta != beta)
	throw(H5::DataSetIException("DensityMatrixPart::load()",
				    "Data in the storage is for another value of the temperature."));
    
    GroundEnergy = Storage->loadReal(RootGroup,"GroundEnergy");
    Z_part = Storage->loadReal(RootGroup,"Z_part");
    Storage->loadRealVector(RootGroup,"weights",weights);
}
