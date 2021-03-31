/** \file include/pomerol/Susceptibility.h
** \brief Dynamical susceptibility.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
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
#include"EnsembleAverage.h"

namespace Pomerol{

/** This class represents a dynamical susceptibility in the Matsubara representation.
 *
 * Exact definition:
 *
 * \f[
 *      \chi(\omega_n) = \int_0^\beta \langle\mathbf{T} A(\tau) B(0)\rangle e^{i\omega_n\tau} d\tau
 * \f]
 *
 * or
 *
 * \f[
 *      \tilde{\chi}(\omega_n) = \chi(\omega_n) - \beta \langle A \rangle \langle B \rangle
 * \f]
 *
 * if specified.
 *
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). A pair of parts, one part of an annihilation operator and
 * another from a creation operator, corresponds to a part of the Green's function.
 */
template<bool Complex = false>
class Susceptibility : public Thermal, public ComputableObject {

public:

    using PartT = SusceptibilityPart<Complex>;

private:

    /** A reference to a states classification object. */
    const StatesClassification<Complex>& S;
    /** A reference to a Hamiltonian. */
    const Hamiltonian<Complex>& H;
    /** A reference to a quadratic operator. */
    const QuadraticOperator<Complex>& A;
    /** A reference to a quadratic operator. */
    const QuadraticOperator<Complex>& B;
    /** A reference to a density matrix. */
    const DensityMatrix<Complex>& DM;

    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool Vanishing;

    /** A list of pointers to parts (every part corresponds to a part of the quadratic operator A
     * and a part of the quadratic operator B).
     */
    std::list<PartT*> parts;

    /** Subtract disconnected part <A><B> */
    bool SubtractDisconnected;

    /** <A>, <B> */
    ComplexType ave_A, ave_B;

public:
     /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] A A reference to a quadratic operator.
     * \param[in] B A reference to a quadratic operator.
     * \param[in] DM A reference to a density matrix.
     */
    Susceptibility(const StatesClassification<Complex>& S,
                   const Hamiltonian<Complex>& H,
                   const QuadraticOperator<Complex>& A,
                   const QuadraticOperator<Complex>& B,
                   const DensityMatrix<Complex>& DM);
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

    /** Activate subtraction of the disconnected part <A><B>
     * <A> and <B> are computed in this class.
     */
    void subtractDisconnected();
    /** Activate subtraction of the disconnected part <A><B>
     * \param[in] ave_A Precomputed value of <A>
     * \param[in] ave_B Precomputed value of <B>
     */
    void subtractDisconnected(ComplexType ave_A, ComplexType ave_B);
    /** Activate subtraction of the disconnected part <A><B>
     * \param[in] EA_A Predefined EnsembleAverage class for operator A.
     * \param[in] EA_B Predefined EnsembleAverage class for operator B.
     */
    void subtractDisconnected(EnsembleAverage<Complex> &EA_A, EnsembleAverage<Complex> &EA_B);

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
template<bool Complex>
inline ComplexType Susceptibility<Complex>::operator()(long int MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing*RealType(2*MatsubaraNumber)); }

template<bool Complex>
inline ComplexType Susceptibility<Complex>::operator()(ComplexType z) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        for(auto iter = parts.begin(); iter != parts.end(); iter++)
            Value += (**iter)(z);
    }
    if(SubtractDisconnected)
        if( abs(z) < 1e-15 )  Value -= ave_A * ave_B * beta;  // only for n=0
    return Value;
}

template<bool Complex>
inline ComplexType Susceptibility<Complex>::of_tau(RealType tau) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        for(auto iter = parts.begin(); iter != parts.end(); iter++)
            Value += (*iter)->of_tau(tau);
    }
    if(SubtractDisconnected)
        Value -= ave_A * ave_B;
    return Value;
}

extern template class Susceptibility<false>;
extern template class Susceptibility<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_SUSCEPTIBILITY_H
