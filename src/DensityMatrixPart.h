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


/** \file src/DensityMatrixPart.h
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_DENSITYMATRIXPART_H
#define __INCLUDE_DENSITYMATRIXPART_H

#include "HDF5Storage.h"
#include "Thermal.h"
#include "HamiltonianPart.h"

namespace Pomerol{

/** This class represents a part of a density matrix.
 * It is responsible for calculations of exponential weights
 * and a contribution to the partition function made by the corresponding
 * block of the Hamiltonian.
 */

class DensityMatrixPart : public HDF5Storable, public Thermal
{
    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a part of a Hamiltonian. */
    const HamiltonianPart& hpart;

    /** The ground energy of the Hamiltonian.
     * It is subtracted from all energy levels to avoid
     * large non-normalized weights and possible loss of precision.
     */
    RealType GroundEnergy;

    /** A real vector holding all weights in this part. */
    RealVectorType weights;

    /** The contribution to the partition function. */
    RealType Z_part;

public:
    /** Constructor.
     * \param[in] hpart A reference to a part of the Hamiltonian.
     * \param[in] beta The inverse temperature.
     * \param[in] GroundEnergy The ground state energy of the Hamiltonian.
     */
    DensityMatrixPart(StatesClassification &S, const HamiltonianPart& hpart, RealType beta, RealType GroundEnergy);

    /** Compute unnormalized weights (diagonal matrix elements)
     * 
     * \return The partial partition function.
     */
    RealType computeUnnormalized(void);

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
    RealType getAverageEnergy(void) const;

    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const;

    /** Returns the partition function of this part. */
    RealType getPartialZ(void) const;

    void save(H5::CommonFG* RootGroup) const;
    void load(const H5::CommonFG* RootGroup);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
