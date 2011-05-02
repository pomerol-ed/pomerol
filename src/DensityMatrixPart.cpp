//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


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
            NN += weights(m)*S.n_i(S.getQuantumState(hpart.id(),fi),i)*S.n_i(S.getQuantumState(hpart.id(),fi),j)*CurrentEigenState(fi)*CurrentEigenState(fi);
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

void DensityMatrixPart::save(H5::CommonFG* RootGroup) const
{
    HDF5Storage::saveReal(RootGroup,"beta",beta);
    HDF5Storage::saveReal(RootGroup,"GroundEnergy",GroundEnergy);
    HDF5Storage::saveReal(RootGroup,"Z_part",Z_part);
    HDF5Storage::saveRealVector(RootGroup,"weights",weights);
}

void DensityMatrixPart::load(const H5::CommonFG* RootGroup)
{
    RealType newBeta = HDF5Storage::loadReal(RootGroup,"beta");
    if(newBeta != beta)
	throw(H5::DataSetIException("DensityMatrixPart::load()",
 				     "Data in the storage is for another value of the temperature."));
    
    GroundEnergy = HDF5Storage::loadReal(RootGroup,"GroundEnergy");
    Z_part = HDF5Storage::loadReal(RootGroup,"Z_part");
    HDF5Storage::loadRealVector(RootGroup,"weights",weights);
}
