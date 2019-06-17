/** \file include/pomerol/Susceptibility.h
** \brief Thermal Green's function.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef __INCLUDE_SUSCEPTIBILITY_H
#define __INCLUDE_SUSCEPTIBILITY_H

#include <sstream>

#include"Misc.h"
#include"Thermal.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"FieldOperator.h"
#include"DensityMatrix.h"
#include"SusceptibilityPart.h"

namespace Pomerol{

/** This class represents a thermal Green's function in the Matsubara representation.
 *
 * Exact definition:
 * 
 * \f[
 *      \chi(\omega_n) = \int_0^\beta \langle\mathbf{T} A(\tau) B(0)\rangle e^{i\omega_n\tau} d\tau
 * \f]
 * 
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). A pair of parts, one part of an annihilation operator and
 * another from a creation operator, corresponds to a part of the Green's function.
 */
class Susceptibility : public Thermal, public ComputableObject {

    /** A reference to a states classification object. */
    const StatesClassification& S;
    /** A reference to a Hamiltonian. */
    const Hamiltonian& H;
    /** A reference to a quadratic operator. */
    const QuadraticOperator& A;
    /** A reference to a quadratic operator. */
    const QuadraticOperator& B;
    /** A reference to a density matrix. */
    const DensityMatrix& DM;

    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool Vanishing;

    /** A list of pointers to parts (every part corresponds to a part of the quadratic operator A
     * and a part of the quadratic operator B).
     */
    std::list<SusceptibilityPart*> parts;

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] A A reference to an operator.
     * \param[in] B A reference to an operator.
     * \param[in] DM A reference to a density matrix.
     */
    Susceptibility(const StatesClassification& S, const Hamiltonian& H,
                   const QuadraticOperator& A, const QuadraticOperator& B, const DensityMatrix& DM);
    /** Copy-constructor.
     * \param[in] GF Susceptibility object to be copied.
     */
    Susceptibility(const Susceptibility& Chi);
    /** Destructor. */
    ~Susceptibility();

    /** Chooses relevant parts of A and B and allocates resources for the parts of the Green's function. */
    void prepare(void);
    /** Actually computes the parts and fills the internal cache of precomputed values.
     * \param[in] NumberOfMatsubaras Number of positive Matsubara frequencies.
     */
    void compute();

    /** Returns the 'bit' (index) of the quadratic operator A or B.
     * \param[in] Position Use A for Position==0 and B for Position==1.
     */
    unsigned short getIndex(size_t Position) const;

     /** Returns the value of the Green's function calculated at a given frequency.
     * \param[in] MatsubaraNum Number of the Matsubara frequency (\f$ \omega_n = \pi(2n+1)/\beta \f$).
     */
    ComplexType operator()(long MatsubaraNumber) const;

     /** Returns the value of the Green's function calculated at a given frequency.
     * \param[in] z Input frequency
     */
    ComplexType operator()(ComplexType z) const;

     /** Returns the value of the Green's function calculated at a given imaginary time point.
     * \param[in] tau Imaginary time point.
     */
    ComplexType of_tau(RealType tau) const;

    bool isVanishing(void) const;
};

// BOSON: bosononic Matsubara frequency
inline ComplexType Susceptibility::operator()(long int MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing*RealType(2*MatsubaraNumber)); }

inline ComplexType Susceptibility::operator()(ComplexType z) const {
    if(Vanishing) return 0;
    else {
        ComplexType Value = 0;
        for(std::list<SusceptibilityPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
            Value += (**iter)(z);
        return Value;
    };
}

inline ComplexType Susceptibility::of_tau(RealType tau) const {
    if(Vanishing) return 0;
    else {
        ComplexType Value = 0;
        for(std::list<SusceptibilityPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
            Value += (*iter)->of_tau(tau);
        return Value;
    };
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_SUSCEPTIBILITY_H

