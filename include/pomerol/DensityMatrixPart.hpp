//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/DensityMatrixPart.hpp
/// \brief Diagonal block of a many-body Gibbs density matrix.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP
#define POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP

#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

namespace Pomerol {

/// \addtogroup ED
///@{

/// \brief Part of a many-body Gibbs density matrix.
///
/// This class represents a part (diagonal block \f$B\f$) of a many-body Gibbs
/// density matrix
/// \f[
///  \hat\rho = \frac{e^{-\beta\hat H}}{Z}, \quad Z = Tr[e^{-\beta\hat H}].
/// \f]
/// Since the matrix is always computed in the eigenbasis of \f$\hat H\f$,
/// it is sufficient to store its eigenvalues (statistical weights
/// \f$w_s = e^{-\beta E_s} / Z\f$).
class DensityMatrixPart : public Thermal {

    /// A reference to the respective diagonal block of the Hamiltonian.
    HamiltonianPart const& H;

    /// The ground state energy of the Hamiltonian.
    /// It is subtracted from all energy levels to avoid exponentially large
    /// unnormalized weights, which could lead to a precision loss during
    /// calculation of the normalized weights.
    RealType GroundEnergy;

    /// Eigenvalues of the density matrix within this block.
    RealVectorType weights;

    /// Contribution of this block to the partition function.
    RealType Z_part = 0;

    /// true, if there are non-negligible weights in this block.
    bool Retained = true;

public:
    /// Constructor.
    /// \param[in] H The respective diagonal block of the Hamiltonian.
    /// \param[in] beta Inverse temperature \f$\beta\f$.
    /// \param[in] GroundEnergy The ground state energy of the Hamiltonian.
    DensityMatrixPart(HamiltonianPart const& H, RealType beta, RealType GroundEnergy);

    /// Compute and store the unnormalized statistical weights \f$Z w_s\f$.
    RealType computeUnnormalized();

    /// Normalize the stored statistical weights by the partition function \f$Z\f$.
    /// \param[in] Z The partition function.
    void normalize(RealType Z);

    /// Return a statistical weight \f$w_s\f$.
    /// \param[in] s Index of the weight within this block.
    RealType getWeight(InnerQuantumState s) const;

    /// Compute the energy averaged over this block, \f$ \langle E\rangle = \sum_{s\in B} E_s w_s\f$.
    RealType getAverageEnergy() const;

    /// Return the contribution of this block to the partition function, \f$Z_B = Z\sum_{s\in B} w_s\f$.
    RealType getPartialZ() const { return Z_part; }

    /// Check if any of the statistical weights is above a given tolerance.
    /// If not, mark this block as irrelevant (set the Retained flag to false).
    /// \param[in] Tolerance Statistical weights smaller or equal to this value are considered negligible.
    void truncate(RealType Tolerance);

    /// Does this block contain any non-negligible statistical weights?
    bool isRetained() const { return Retained; }
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP
