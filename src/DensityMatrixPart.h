/** \file src/DensityMatrixPart.h
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_DENSITYMATRIXPART_H
#define __INCLUDE_DENSITYMATRIXPART_H

#include "HDF5Storage.h"
#include "HamiltonianPart.h"

/** This class represents a part of a density matrix.
 * It is responsible for calculations of exponential weights
 * and a contribution to the partition function made by the corresponding
 * block of the Hamiltonian.
 */

class DensityMatrixPart
#ifdef pomerolHDF5
: public HDF5Storable
#endif
{
    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a part of a Hamiltonian. */
    HamiltonianPart& hpart;
    /** The inverse temperature. */
    RealType beta;
    /** The ground energy of the Hamiltonian.
     * It is subtracted from all energy levels to avoid
     * large non-normalized weights and possible loss of precision.
     */
    RealType GroundEnergy;

    /** The number of states in this part. */
    QuantumState partSize;
    /** A real vector holding all weights in this part. */
    RealVectorType weights;

    /** A contribution to the partition function. */
    RealType Z_part;

public:
    /** Constructor.
     * \param[in] hpart A reference to a part of the Hamiltonian.
     * \param[in] beta The inverse temperature.
     * \param[in] GroundEnergy The ground state energy of the Hamiltonian.
     */
    DensityMatrixPart(StatesClassification &S, HamiltonianPart& hpart, RealType beta, RealType GroundEnergy);
    /** Divide all the weights by the partition function.
     * 
     * Warning! Must be called once and only once!
     * \param[in] Z The partition function.
     */
    void normalize(RealType Z);

    /** Returns an averaged value of the energy. */
    RealType getAverageEnergy(void);
    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j);
    /** Performs computations of the weights and a contribution to the partition function. */
    RealType compute(void);
    /** Returns the weight corresponding to a specified number of state.
     * \param[in] m A number of a state inside this part.
     */
    RealType weight(int m);
    /** Returns the inverse temperature. */
    RealType getBeta();

    /** Returns the number of the state, which has the lowest statistical weight before exceeding TruncationTolerance 
     * \param[in] TruncationTolerance - the level at which the statistical weight should be cutted
     */
    InnerQuantumState getMaximumTruncationState( RealType TruncationTolerance);
#ifdef pomerolHDF5 
    void save(H5::CommonFG* FG) const;
    void load(const H5::CommonFG* FG);
#endif
};

#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
