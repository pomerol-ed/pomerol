/** \file src/DensityMatrix.h
** \brief Density matrix of the grand canonical ensemble.
** 
** \author Igor Krivenko (igor@shg.ru)
*/

#ifndef ____DEFINE_DENSITY_MATRIX____
#define ____DEFINE_DENSITY_MATRIX____

#include "StatesClassification.h"
#include "BitClassification.h"
#include "Hamiltonian.h"
#include "DensityMatrixPart.h"

/** This class represents a density matrix /f$\rho = \exp(-\beta \hat H)/Z$/f.
 * It is actually a container class for a collection of parts (all real calculations
 * take place inside the parts). There is one-to-one correspondence between parts of
 * the Hamiltonian and the parts of the density matrix itself, since the density matrix
 * is a function of /f\hat H/f.
 */
class DensityMatrix
{
    /** Reference to a states classification object. */
    StatesClassification& S;
    /** Reference to a Hamiltonian defining the grand canonical ensemble. */
    Hamiltonian &H;
    /** Array of pointers to parts (every part corresponds to a part of the Hamiltonian). */
    DensityMatrixPart** parts;
    /** Inverse temperature. */
    RealType beta;

    /** Number of parts. */
    BlockNumber NumOfBlocks;

public:
    /** Constructor.
     * \param[in] S Reference to a states classification object.
     * \param[in] H Reference to a Hamiltonian.
     * \param[in] beta Inverse temperature.
     */
    DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta);
    /** Destructor. */
    ~DensityMatrix();

    /** Returns a part of the density matrix.
    * \param[in] in Set of the quantum numbers to be resolved into a part number.
    */
    DensityMatrixPart& part(const QuantumNumbers &in);
    /** Returns a part of the density matrix.
     * \param[in] in Part number.
     */
    DensityMatrixPart& part(BlockNumber in);

    /** Allocates resources for the parts. */
    void prepare(void);
    /** Actually computes the parts. */
    void compute(void);
    /** Returns the value of the density matrix corresponding to a specified quantum state.
     * \param[in] state Quantum state.
     */
    RealType operator()(QuantumState &state);
    /** Returns the inverse temperature. */
    RealType getBeta();
    /** Returns an averaged value of the energy. */
	RealType getAverageEnergy();
};

#endif // endif :: #ifndef ____DEFINE_DENSITY_MATRIX____
