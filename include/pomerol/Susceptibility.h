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
 *      G(\omega_n) = -\int_0^\beta \langle\mathbf{T}c_i(\tau)c^+_j(0)\rangle e^{i\omega_n\tau} d\tau
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
    /** A reference to an annihilation operator. */
    const AnnihilationOperator& C;
    /** A reference to a creation operator. */
    const CreationOperator& CX;
    /** A reference to a density matrix. */
    const DensityMatrix& DM;

    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool Vanishing;

    /** A list of pointers to parts (every part corresponds to a part of the annihilation operator
     * and a part of the creation operator).
     */
    std::list<SusceptibilityPart*> parts;

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] C A reference to an annihilation operator.
     * \param[in] CX A reference to a creation operator.
     * \param[in] DM A reference to a density matrix.
     */
    Susceptibility(const StatesClassification& S, const Hamiltonian& H,
                   const AnnihilationOperator& C, const CreationOperator& CX, const DensityMatrix& DM);
    /** Copy-constructor.
     * \param[in] GF Susceptibility object to be copied.
     */
    Susceptibility(const Susceptibility& GF);
    /** Destructor. */
    ~Susceptibility();

    /** Chooses relevant parts of C and CX and allocates resources for the parts of the Green's function. */
    void prepare(void);
    /** Actually computes the parts and fills the internal cache of precomputed values.
     * \param[in] NumberOfMatsubaras Number of positive Matsubara frequencies.
     */
    void compute();

    /** Returns the 'bit' (index) of the operator C or CX.
     * \param[in] Position Use C for Position==0 and CX for Position==1.
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

inline ComplexType Susceptibility::operator()(long int MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing*RealType(2*MatsubaraNumber+1)); }

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

