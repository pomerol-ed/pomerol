/** \file include/pomerol/EnsembleAverage.h
** \brief Thermal Green's function.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_ENSEMBLEAVERAGE_H
#define __INCLUDE_ENSEMBLEAVERAGE_H

#include <sstream>

#include"Misc.h"
#include"Thermal.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"FieldOperator.h"
#include"DensityMatrix.h"
#include"EnsembleAveragePart.h"

namespace Pomerol{

/** This class represents the ensemble average of a bosonic operator.
 *
 * Exact definition:
 * 
 * \f[
 *      \langle A \rangle
 * \f]
 * 
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). A pair of parts, one part of an annihilation operator and
 * another from a creation operator, corresponds to a part of the Green's function.
 */
class EnsembleAverage : public Thermal, public ComputableObject {

    /** A reference to a states classification object. */
    const StatesClassification& S;
    /** A reference to a Hamiltonian. */
    const Hamiltonian& H;
    /** A reference to a bosonic operator. */
    const QuadraticOperator& A;
    /** A reference to a density matrix. */
    const DensityMatrix& DM;

    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool Vanishing;

    /** A list of pointers to parts (every part corresponds to a part of the annihilation operator
     * and a part of the creation operator).
     */
    std::list<EnsembleAveragePart*> parts;

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] A A reference to a quadratic operator.
     * \param[in] DM A reference to a density matrix.
     */
     EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                     const QuadraticOperator& A, const DensityMatrix& DM);
    /** Copy-constructor.
     * \param[in] EA EnsembleAverage object to be copied.
     */
    EnsembleAverage(const EnsembleAverage& EA);
    /** Destructor. */
    ~EnsembleAverage();

    /** Chooses relevant parts of A and allocates resources for the parts of the ensemble average. */
    void prepare(void);
    /** Actually computes the parts and fills the internal cache of precomputed values. */
    void compute();

    ComplexType getResult(){
        if(Vanishing) return 0;
        else {
            ComplexType Value = 0;
            for(std::list<EnsembleAveragePart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
                Value += (**iter).getResult();
            return Value;
        };
    };

    bool isVanishing(void) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_ENSEMBLEAVERAGE_H
