//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/DensityMatrixPart.h
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP
#define POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP

#include "HamiltonianPart.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"

namespace Pomerol {

/** This class represents a part of a density matrix.
 * It is responsible for calculations of exponential weights
 * and a contribution to the partition function made by the corresponding
 * block of the Hamiltonian.
 */

class DensityMatrixPart : public Thermal {
    /** A reference to a part of a Hamiltonian. */
    HamiltonianPart const& H;

    /** The ground energy of the Hamiltonian.
     * It is subtracted from all energy levels to avoid
     * large non-normalized weights and possible loss of precision.
     */
    RealType GroundEnergy;

    /** A real vector holding all weights in this part. */
    RealVectorType weights;

    /** The contribution to the partition function. */
    RealType Z_part = 0;

    /** It is true if this part has not been truncated. */
    bool Retained = true;

public:
    /** Constructor.
     * \param[in] hpart A reference to a part of the Hamiltonian.
     * \param[in] beta The inverse temperature.
     * \param[in] GroundEnergy The ground state energy of the Hamiltonian.
     */
    DensityMatrixPart(HamiltonianPart const& H, RealType beta, RealType GroundEnergy);

    /** Compute unnormalized weights (diagonal matrix elements)
     *
     * \return The partial partition function.
     */
    RealType computeUnnormalized();

    /** Divide all the weights by the partition function.
     *
     * \param[in] Z The partition function.
     */
    void normalize(RealType Z);

    /** Returns the weight corresponding to a specified state.
     * \param[in] s State inside this part.
     */
    RealType getWeight(InnerQuantumState s) const;

    /** Returns an averaged value of the energy. */
    RealType getAverageEnergy() const;

    /** Returns the partition function of this part. */
    RealType getPartialZ() const { return Z_part; }

    /** Truncates this part if it does not include any states having larger weight than Tolerance. */
    void truncate(RealType Tolerance);

    /** Returns true if this part has not been truncated. */
    bool isRetained() const { return Retained; }
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_DENSITYMATRIXPART_HPP
