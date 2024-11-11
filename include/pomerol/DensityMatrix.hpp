//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/DensityMatrix.hpp
/// \brief Many-body Gibbs density matrix as a list of diagonal blocks.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_HPP
#define POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_HPP

#include "ComputableObject.hpp"
#include "DensityMatrixPart.hpp"
#include "Hamiltonian.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <vector>

namespace Pomerol {

/// \addtogroup ED
///@{

/// \brief Many-body Gibbs density matrix.
///
/// This class represents a many-body Gibbs density matrix
/// \f[
///  \hat\rho = \frac{e^{-\beta\hat H}}{Z}, \quad Z = Tr[e^{-\beta\hat H}].
/// \f]
/// The matrix is stored as a list of \ref DensityMatrixPart (diagonal blocks), which correspond
/// to invariant subspaces/diagonal blocks of the Hamiltonian \f$\hat H\f$.
class DensityMatrix : public Thermal, public ComputableObject {

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// A reference to the Hamiltonian \f$\hat H\f$.
    Hamiltonian const& H;
    /// The list of parts (diagonal blocks).
    std::vector<DensityMatrixPart> parts;

public:
    /// Constructor.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] H The Hamiltonian \f$\hat H\f$.
    /// \param[in] beta Inverse temperature \f$\beta\f$.
    DensityMatrix(StatesClassification const& S, Hamiltonian const& H, RealType beta);

    /// Destructor.
    ~DensityMatrix() = default;

    /// Allocate memory for the parts.
    void prepare();

    /// Compute statistical weights within every part (diagonal block).
    /// \pre prepare() has been called.
    void compute();

    /// Return a reference to a part (diagonal block).
    /// \param[in] B Index of the part.
    DensityMatrixPart const& getPart(BlockNumber B) const;

    /// Return a statistical weight corresponding to a specified eigenstate.
    /// \param[in] state Index of the eigenstate within the full Hilbert space.
    RealType getWeight(QuantumState state) const;

    /// Compute the average energy \f$ \langle E\rangle = \sum_s E_s w_s\f$.
    RealType getAverageEnergy() const;

    /// Check if any of the statistical weights within each block is above a given tolerance.
    /// If not, mark the respective block as irrelevant (set the Retained flag to false).
    /// \param[in] Tolerance Statistical weights smaller or equal to this value are considered negligible.
    /// \param[in] verbose Print out information about the truncation results.
    void truncateBlocks(RealType Tolerance, bool verbose = true);

    /// Does a given block contain any non-negligible statistical weights?
    /// \param[in] B Index of the part (block).
    bool isRetained(BlockNumber B) const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_HPP
