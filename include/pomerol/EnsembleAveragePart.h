/** \file include/pomerol/EnsembleAveragePart.h
** \brief Part of a Green's function for a given set of quantum numbers.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_ENSEMBLEAVERAGEPART_H
#define __INCLUDE_ENSEMBLEAVERAGEPART_H

#include<iomanip>
#include<cmath>

#include"Misc.h"
#include"StatesClassification.h"
#include"HamiltonianPart.h"
#include"FieldOperator.h"
#include"DensityMatrixPart.h"
#include"TermList.h"

namespace Pomerol{

/** This class represents a part of an ensemble average.
 * Every part describes all transitions allowed by selection rules
 * between a given pair of Hamiltonian blocks.
 */
class EnsembleAveragePart : public Thermal
{
    /** A reference to a part of a Hamiltonian. */
    const HamiltonianPart& Hpart;
    /** A reference to a part of a density matrix. */
    const DensityMatrixPart& DMpart;

    /** A reference to a part of an annihilation operator. */
    const QuadraticOperatorPart& A;

    /** The result is summed up here **/
    ComplexType result;

    /** A matrix element with magnitude less than this value is treated as zero. */
    const RealType MatrixElementTolerance; // 1e-8;

public:

    /** Constructor.
     * \param[in] A A reference to a part of an annihilation operator.
     * \param[in] Hpart A reference to a part of the Hamiltonian.
     * \param[in] DMpart A reference to a part of the density matrix.
     */
    EnsembleAveragePart(const QuadraticOperatorPart& A,
                       const HamiltonianPart& Hpart, const DensityMatrixPart& DMpart);

    /** Iterates over all matrix elements and fills the list of terms. */
    void compute(void);

    /** Get result for ensemble average */
    ComplexType getResult(){ return result;};
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_ENSEMBLEAVERAGEPART_H
