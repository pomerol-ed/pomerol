//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#include "DensityMatrixPart.h"

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
