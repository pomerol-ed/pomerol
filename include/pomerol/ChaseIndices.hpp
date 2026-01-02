//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/ChaseIndices.hpp
/// \brief `Chase indices` algorithm used in multipoint correlator calculations.
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_POMEROL_CHASEINDICES_HPP
#define POMEROL_INCLUDE_POMEROL_CHASEINDICES_HPP

#include "pomerol/StatesClassification.hpp"

#ifndef DOXYGEN_SKIP
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Pomerol {

/// Make the lagging index catch up or outrun the leading index.
template <bool Complex>
inline bool chaseIndices(typename RowMajorMatrixType<Complex>::InnerIterator& index1_iter,
                         typename ColMajorMatrixType<Complex>::InnerIterator& index2_iter) {
    InnerQuantumState index1 = index1_iter.index();
    InnerQuantumState index2 = index2_iter.index();

    if(index1 == index2)
        return true;

    if(index1 < index2)
        for(; InnerQuantumState(index1_iter.index()) < index2 && index1_iter; ++index1_iter)
            ;
    else
        for(; InnerQuantumState(index2_iter.index()) < index1 && index2_iter; ++index2_iter)
            ;

    return false;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_CHASEINDICES_HPP
