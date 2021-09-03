//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/DensityMatrix.h
** \brief Density matrix of the grand canonical ensemble.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_H
#define POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_H

#include "ComputableObject.hpp"
#include "DensityMatrixPart.hpp"
#include "Hamiltonian.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

#include <vector>

namespace Pomerol {

/** This class represents a density matrix \f$ \rho = \exp(-\beta \hat H)/Z \f$.
 * It is actually a container class for a collection of parts (all real calculations
 * take place inside the parts). There is one-to-one correspondence between parts of
 * the Hamiltonian and the parts of the density matrix itself, since the density matrix
 * is a function of \f$ \hat H \f$.
 */
class DensityMatrix : public Thermal, public ComputableObject {
    /** A reference to a states classification object. */
    StatesClassification const& S;
    /** A reference to a Hamiltonian defining the grand canonical ensemble. */
    Hamiltonian const& H;
    /** A vector of pointers to parts (every part corresponds to a part of the Hamiltonian). */
    std::vector<DensityMatrixPart> parts;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] beta The inverse temperature.
     */
    DensityMatrix(StatesClassification const& S, Hamiltonian const& H, RealType beta);
    /** Destructor. */
    ~DensityMatrix() = default;

    /** Allocates resources for the parts. */
    void prepare();

    /** Actually computes the parts. */
    void compute();

    /** Returns a part of the density matrix.
     * \param[in] in A part number.
     */
    DensityMatrixPart const& getPart(BlockNumber in) const;

    /** Returns the value of the density matrix corresponding to a specified quantum state.
     * \param[in] state A quantum state.
     */
    RealType getWeight(QuantumState state) const;

    /** Returns the average energy. */
    RealType getAverageEnergy() const;

    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const;

    /** Truncate such blocks that do not include any states having larger weight than Tolerance. */
    void truncateBlocks(RealType Tolerance, bool verbose = true);

    /** Return true if the block has not been truncated. Always true if function truncateBlocks has not been called. */
    bool isRetained(BlockNumber in) const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIX_H
