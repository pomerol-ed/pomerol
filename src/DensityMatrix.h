// This file is part of pomerol ED code
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


/** \file src/DensityMatrix.h
** \brief Density matrix of the grand canonical ensemble.
** 
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_DENSITYMATRIX_H
#define __INCLUDE_DENSITYMATRIX_H

#include "HDF5Storage.h"
#include "ComputableObject.h"
#include "StatesClassification.h"
#include "IndexClassification.h"
#include "Hamiltonian.h"
#include "DensityMatrixPart.h"

/** This class represents a density matrix \f$ \rho = \exp(-\beta \hat H)/Z \f$.
 * It is actually a container class for a collection of parts (all real calculations
 * take place inside the parts). There is one-to-one correspondence between parts of
 * the Hamiltonian and the parts of the density matrix itself, since the density matrix
 * is a function of \f$ \hat H \f$.
 */

class DensityMatrix : public ComputableObject, public Thermal, public HDF5Storable
{
    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a Hamiltonian defining the grand canonical ensemble. */
    Hamiltonian &H;
    /** A vector of pointers to parts (every part corresponds to a part of the Hamiltonian). */
    std::vector<DensityMatrixPart*> parts;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] beta The inverse temperature.
     */
    DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta);
    /** Destructor. */
    ~DensityMatrix();

    /** Returns a part of the density matrix.
    * \param[in] in A set of the quantum numbers to be resolved into a part number.
    */
    DensityMatrixPart& part(const QuantumNumbers &in);
    /** Returns a part of the density matrix.
     * \param[in] in A part number.
     */
    DensityMatrixPart& part(BlockNumber in);

    /** Allocates resources for the parts. */
    void prepare(void);
    /** Actually computes the parts. */
    void compute(void);
    /** Returns the value of the density matrix corresponding to a specified quantum state.
     * \param[in] state A quantum state.
     */
    RealType operator()(QuantumState &state);
    /** Returns the inverse temperature. */
    RealType getAverageEnergy();
    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j);
    
    void save(H5::CommonFG* RootGroup) const;
    void load(const H5::CommonFG* RootGroup);
};

#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
