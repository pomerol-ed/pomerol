//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/Susceptibility.h
** \brief Dynamical susceptibility.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
*/
#ifndef POMEROL_INCLUDE_SUSCEPTIBILITY_HPP
#define POMEROL_INCLUDE_SUSCEPTIBILITY_HPP

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "EnsembleAverage.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "SusceptibilityPart.hpp"
#include "Thermal.hpp"

#include <complex>
#include <numeric>
#include <vector>

namespace Pomerol {

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
class Susceptibility : public Thermal, public ComputableObject {

    /** A reference to a states classification object. */
    StatesClassification const& S;
    /** A reference to a Hamiltonian. */
    Hamiltonian const& H;
    /** A reference to a quadratic operator. */
    MonomialOperator const& A;
    /** A reference to a quadratic operator. */
    MonomialOperator const& B;
    /** A reference to a density matrix. */
    DensityMatrix const& DM;

    /** A flag to represent if Greens function vanishes, i.e. identical to 0 */
    bool Vanishing = true;

    /** A list of pointers to parts (every part corresponds to a part of the quadratic operator A
     * and a part of the quadratic operator B).
     */
    std::vector<SusceptibilityPart> parts;

    /** Subtract disconnected part <A><B> */
    bool SubtractDisconnected = false;

    /** <A>, <B> */
    ComplexType ave_A = {}, ave_B = {};

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] A A reference to a quadratic operator.
     * \param[in] B A reference to a quadratic operator.
     * \param[in] DM A reference to a density matrix.
     */
    Susceptibility(StatesClassification const& S,
                   Hamiltonian const& H,
                   MonomialOperator const& A,
                   MonomialOperator const& B,
                   DensityMatrix const& DM);
    /** Copy-constructor.
     * \param[in] GF Susceptibility object to be copied.
     */
    Susceptibility(Susceptibility const& Chi);

    /** Chooses relevant parts of A and B and allocates resources for the parts of the Green's function. */
    void prepare();
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
    void subtractDisconnected(EnsembleAverage& EA_A, EnsembleAverage& EA_B);

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

    bool isVanishing() const { return Vanishing; }
};

// BOSON: bosononic Matsubara frequency
inline ComplexType Susceptibility::operator()(long int MatsubaraNumber) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber));
}

inline ComplexType Susceptibility::operator()(ComplexType z) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        Value = std::accumulate(parts.begin(),
                                parts.end(),
                                ComplexType(0),
                                [z](ComplexType s, SusceptibilityPart const& p) { return s + p(z); });
    }
    if(SubtractDisconnected && std::abs(z) < 1e-15)
        Value -= ave_A * ave_B * beta; // only for n=0
    return Value;
}

inline ComplexType Susceptibility::of_tau(RealType tau) const {
    ComplexType Value = 0;
    if(!Vanishing) {
        Value = std::accumulate(parts.begin(),
                                parts.end(),
                                ComplexType(0),
                                [tau](ComplexType s, SusceptibilityPart const& p) { return s + p.of_tau(tau); });
    }
    if(SubtractDisconnected)
        Value -= ave_A * ave_B;
    return Value;
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_SUSCEPTIBILITY_HPP
