/** \file src/DensityMatrix.h
** \brief Density matrix of the grand canonical ensemble.
** 
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_DENSITYMATRIX_H
#define __INCLUDE_DENSITYMATRIX_H

#include "Dumper.h"
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
#ifdef pomerolHDF5
class DensityMatrix : public Dumpable
#else
class DensityMatrix
#endif
{
    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a Hamiltonian defining the grand canonical ensemble. */
    Hamiltonian &H;
    /** A vector of pointers to parts (every part corresponds to a part of the Hamiltonian). */
    std::vector<DensityMatrixPart*> parts;
    /** The inverse temperature. */
    RealType beta;

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
    RealType getBeta() const;
    /** Returns an averaged value of the energy. */
    RealType getAverageEnergy();
    
#ifdef pomerolHDF5
    void dumpIt(H5::CommonFG* FG) const;
#endif
};

#endif // endif :: #ifndef __INCLUDE_DENSITYMATRIXPART_H
