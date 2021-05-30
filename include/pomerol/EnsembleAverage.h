/** \file include/pomerol/EnsembleAverage.h
** \brief Ensemble average.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
*/
#ifndef POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_H
#define POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_H

#include"Misc.h"
#include"Thermal.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"MonomialOperator.h"
#include"DensityMatrix.h"

namespace Pomerol {

/** This class represents the ensemble average of a quadratic operator.
 *
 * Exact definition:
 *
 * \f[
 *      \langle A \rangle = \langle c_i^{\dag} c_j \rangle
 * \f]
 *
 * How to use:
 *   EnsembleAverage EA(S, H, A, DM);
 *   EA.prepare()
 *   EA.getResult()
 */
class EnsembleAverage : public Thermal, public ComputableObject {

    /** A reference to a states classification object. */
    const StatesClassification& S;
    /** A reference to a Hamiltonian. */
    const Hamiltonian& H;
    /** A reference to a bosonic operator. */
    const MonomialOperator& A;
    /** A reference to a density matrix. */
    const DensityMatrix& DM;

    ComplexType Result = 0;

    /** Returns the contribution to the ensemble average from a part. Called in prepare() */
    template<bool Complex>
    ComplexType computeImpl(const MonomialOperatorPart& Apart, const DensityMatrixPart& DMpart);

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] A A reference to a quadratic operator.
     * \param[in] DM A reference to a density matrix.
     */
     EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                     const MonomialOperator& A, const DensityMatrix& DM);
    /** Copy-constructor.
     * \param[in] EA EnsembleAverage object to be copied.
     */
    EnsembleAverage(const EnsembleAverage& EA);

    /** Compute the ensemble average of A by choosing relevant parts of A and sum up each contribution. */
    void compute();

    /** Returns the ensemble average */
    ComplexType getResult() const { return Result; };
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_ENSEMBLEAVERAGE_H
