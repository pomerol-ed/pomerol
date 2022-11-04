//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/DensityMatrixPart.cpp
/// \brief Diagonal block of a many-body Gibbs density matrix (implementation).
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/DensityMatrixPart.hpp"

namespace Pomerol {

DensityMatrixPart::DensityMatrixPart(HamiltonianPart const& H, RealType beta, RealType GroundEnergy)
    : Thermal(beta), H(H), GroundEnergy(GroundEnergy), weights(H.getSize()) {}

RealType DensityMatrixPart::computeUnnormalized() {
    weights = exp(-beta * (H.getEigenValues().array() - GroundEnergy));
    return weights.sum();
}

void DensityMatrixPart::normalize(RealType Z) {
    weights /= Z;
    Z_part /= Z;
}

RealType DensityMatrixPart::getAverageEnergy() const {
    return weights.dot(H.getEigenValues());
}

RealType DensityMatrixPart::getWeight(InnerQuantumState s) const {
    return weights(static_cast<Eigen::Index>(s));
}

void DensityMatrixPart::truncate(RealType Tolerance) {
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
