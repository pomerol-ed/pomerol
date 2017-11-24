/** \file include/pomerol/DensityMatrix.h
** \brief Density matrix of the grand canonical ensemble.
** 
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_DENSITYMATRIX_H
#define __INCLUDE_DENSITYMATRIX_H

#include "StatesClassification.h"
#include "IndexClassification.h"
#include "Hamiltonian.h"
#include "DensityMatrixPart.h"

namespace Pomerol{

/** This class represents a density matrix \f$ \rho = \exp(-\beta \hat H)/Z \f$.
 * It is actually a container class for a collection of parts (all real calculations
 * take place inside the parts). There is one-to-one correspondence between parts of
 * the Hamiltonian and the parts of the density matrix itself, since the density matrix
 * is a function of \f$ \hat H \f$.
 */
class DensityMatrix : public Thermal, public ComputableObject
{
    /** A reference to a states classification object. */
    const StatesClassification& S;
    /** A reference to a Hamiltonian defining the grand canonical ensemble. */
    const Hamiltonian &H;
    /** A vector of pointers to parts (every part corresponds to a part of the Hamiltonian). */
    std::vector<DensityMatrixPart*> parts;
    /** A vector of bool: true if the block has not been truncated. */
    std::vector<bool> block_retained;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] beta The inverse temperature.
     */
    DensityMatrix(const StatesClassification& S, const Hamiltonian& H, RealType beta);
    /** Destructor. */
    ~DensityMatrix();

    /** Allocates resources for the parts. */
    void prepare(void);

    /** Actually computes the parts. */
    void compute(void);

    /** Returns a part of the density matrix.
    * \param[in] in A set of the quantum numbers to be resolved into a part number.
    */
    const DensityMatrixPart& getPart(const QuantumNumbers &in) const;
    /** Returns a part of the density matrix.
     * \param[in] in A part number.
     */
    const DensityMatrixPart& getPart(BlockNumber in) const;

    /** Returns the value of the density matrix corresponding to a specified quantum state.
     * \param[in] state A quantum state.
     */
    RealType getWeight(QuantumState state) const;

    /** Returns the average energy. */
    RealType getAverageEnergy() const;
    /** Returns the total average occupancy. */
    RealType getAverageOccupancy() const;
    /** Returns the average occupancy at site i. */
    RealType getAverageOccupancy(ParticleIndex i) const;

    /** Returns an averaged value of the double occupancy. */
    RealType getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const;

    /** Truncate such blocks that do not include any states having larger weight than Tolerance. */
    void truncateBlocks(RealType Tolerance, bool verbose=true);

    /** Return true if the block has not been truncated. Always true if function truncateBlocks has not been called. */
    bool isRetained(BlockNumber in) const;
};

}; // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
