//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/EnsembleAverage.cpp
/// \brief Ensemble average of a monomial operator representing a physical observable (implementation).
/// \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/EnsembleAverage.hpp"

namespace Pomerol {

EnsembleAverage::EnsembleAverage(MonomialOperator const& A, DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), A(A), DM(DM) {}

EnsembleAverage::EnsembleAverage(EnsembleAverage const& EA)
    : Thermal(EA.beta), ComputableObject(EA), A(EA.A), DM(EA.DM), Result(EA.Result) {}

void EnsembleAverage::compute() {
    if(getStatus() >= Prepared)
        return;

    // Find out non-trivial blocks of A.
    MonomialOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(auto Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++) {
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright) {
            // check if retained blocks are included. If not, do not push.
            if(DM.isRetained(Aleft)) {
                if(A.isComplex())
                    Result += computeImpl<true>(const_cast<MonomialOperatorPart&>(A.getPartFromLeftIndex(Aleft)),
                                                DM.getPart(Aleft));
                else
                    Result += computeImpl<false>(const_cast<MonomialOperatorPart&>(A.getPartFromLeftIndex(Aleft)),
                                                 DM.getPart(Aleft));
            }
        }
    }

    setStatus(Prepared);
}

// This function is called directly in prepare()
template <bool Complex>
ComplexType EnsembleAverage::computeImpl(MonomialOperatorPart const& Apart, DensityMatrixPart const& DMpart) {
    // Blocks (submatrices) of A
    RowMajorMatrixType<Complex> const& Amatrix = Apart.getRowMajorValue<Complex>();

    // Sum up <index1|A|index1> * weight(index1)
    ComplexType result_part = 0;
    for(Eigen::Index Index = 0; Index < Amatrix.outerSize(); ++Index)
        result_part += Amatrix.coeff(Index, Index) * DMpart.getWeight(Index);
    return result_part;
}

} // namespace Pomerol
