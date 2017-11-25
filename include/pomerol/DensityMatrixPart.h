/** \file include/pomerol/DensityMatrixPart.h
** \brief Part (diagonal block) of a density matrix.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_DENSITYMATRIXPART_H
#define __INCLUDE_DENSITYMATRIXPART_H

#include "Thermal.h"
#include "HamiltonianPart.h"

namespace Pomerol{

/** This class represents a part of a density matrix.
 * It is responsible for calculations of exponential weights
 * and a contribution to the partition function made by the corresponding
 * block of the Hamiltonian.
 */

class DensityMatrixPart : public Thermal
{
    /** A reference to a states classification object. */
    const StatesClassification& S;
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

    /** It is true if this part has not been truncated. */
    bool retained;

public:
    /** Constructor.
     * \param[in] hpart A reference to a part of the Hamiltonian.
     * \param[in] beta The inverse temperature.
     * \param[in] GroundEnergy The ground state energy of the Hamiltonian.
     */
    DensityMatrixPart(const StatesClassification &S, const HamiltonianPart& hpart, RealType beta, RealType GroundEnergy);

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

    RealType getAverageOccupancy() const;
    /** Returns the average occupancy at site i. */
    RealType getAverageOccupancy(ParticleIndex i) const;
    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const;

    /** Returns the partition function of this part. */
    RealType getPartialZ(void) const;

    /** Truncates this part if it does not include any states having larger weight than Tolerance. */
    void truncate(RealType Tolerance);

    /** Returns true if this part has not been truncated. */
    bool isRetained() const;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
