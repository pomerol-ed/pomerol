// This file is part of pomerol ED code
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


/** \file src/TwoParticleGFPart.h
** \brief Part of a two-particle Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_TWOPARTICLEGFPART_H
#define __INCLUDE_TWOPARTICLEGFPART_H

#include"Misc.h"
#include"ComputableObject.h"
#include"StatesClassification.h"
#include"HamiltonianPart.h"
#include"FieldOperatorPart.h"
#include"DensityMatrixPart.h"

/** This class represents a part of a two-particle Green's function.
 * Every part describes one 'world stripe' of four operators.
 */
class TwoParticleGFPart : public ComputableObject, public Thermal {

friend class TwoParticleGF;

public:

    /** A non-resonant term has the following form:
     * \f[
     * \frac{C}{(z_1-P_1)(z_2-P_2)(z_3-P_3)}
     * \f]
     * if isz4 == false, and
     * \f[
     * \frac{C}{(z_1-P_1)(z_1+z_2+z_3-P_1-P_2-P_3)(z_3-P_3)}
     * \f]
     * otherwise.
     */
    struct NonResonantTerm;
     
    /** A resonant term has the following form:
    * \f[
    * \frac{1}{(z_1-P_1)(z_3-P_3)}
    *   \left( R \delta(z_1+z_2-P_1-P_2) + N \frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2} \right)
    * \f]
    */
    struct ResonantTerm;
   
private:

    /** A reference to a part of the first operator. */
    FieldOperatorPart& O1;
    /** A reference to a part of the second operator. */
    FieldOperatorPart& O2;
    /** A reference to a part of the third operator. */
    FieldOperatorPart& O3;
    /** A reference to a part of the fourth (creation) operator. */
    CreationOperatorPart& CX4;

    /** A reference to the first part of a Hamiltonian. */
    HamiltonianPart& Hpart1;
    /** A reference to the second part of a Hamiltonian. */
    HamiltonianPart& Hpart2;
    /** A reference to the third part of a Hamiltonian. */
    HamiltonianPart& Hpart3;
    /** A reference to the fourth part of a Hamiltonian. */
    HamiltonianPart& Hpart4;

    /** A reference to the first part of a density matrix (the part corresponding to Hpart1). */
    DensityMatrixPart& DMpart1;
    /** A reference to the second part of a density matrix (the part corresponding to Hpart2). */
    DensityMatrixPart& DMpart2;
    /** A reference to the third part of a density matrix (the part corresponding to Hpart3). */
    DensityMatrixPart& DMpart3;
    /** A reference to the fourth part of a density matrix (the part corresponding to Hpart4). */
    DensityMatrixPart& DMpart4;

    /** A permutation of the operators for this part. */
    Permutation3 Permutation;

    /** A list of non-resonant terms. */
    std::list<NonResonantTerm> NonResonantTerms;
    /** A list of resonant terms. */
    std::list<ResonantTerm> ResonantTerms;
    
    
    /** Adds a multi-term that has the following form:
    * \f[
    * \frac{1}{(z_1-P_1)(z_3-P_3)}
    *         \left(\frac{C_4}{z_1+z_2+z_3-P_1-P_2-P_3} + \frac{C_2}{z_2-P_2} \right. +
    * \f]
    * \f[     \left.
    *         + R_{12}\delta(z_1+z_2-P_1-P_2)
    *         + N_{12}\frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2}
    *         + R_{23}\delta(z_2+z_3-P_2-P_3)
    *         + N_{23}\frac{1 - \delta(z_2+z_3-P_2-P_3)}{z_2+z_3-P_2-P_3}
    *         \right)
    * \f]
    *
    * Where
    * \f{eqnarray*}{
    *      P_1 = E_j - E_i \\
    *      P_2 = E_k - E_j \\
    *      P_3 = E_l - E_k \\
    *      C_2 = -C(w_j + w_k) \\
    *      C_4 = C(w_i + w_l) \\
    *      R_{12} = C\beta w_i \\
    *      N_{12} = C(w_k - w_i) \\
    *      R_{23} = -C\beta w_j \\
    *      N_{23} = C(w_j - w_l)
    * \f}
    *
    * In fact this is a slightly rewritten form of an equation for \f$ \phi \f$ from
    * <em>H. Hafermann et al 2009 EPL 85 27007</em>.
    * 
    * \param[in] Coeff Common prefactor \f$ C \f$ for coefficients \f$ C_2 \f$, \f$ C_4 \f$,
    *              \f$ R_{12} \f$, \f$ N_{12} \f$, \f$ R_{23} \f$, \f$ N_{23} \f$.
    * \param[in] beta The inverse temperature.
    * \param[in] Ei The first energy level \f$ E_i \f$.
    * \param[in] Ej The second energy level \f$ E_j \f$.
    * \param[in] Ek The third energy level \f$ E_k \f$.
    * \param[in] El The fourth energy level \f$ E_l \f$.
    * \param[in] Wi The first weight \f$ w_i \f$.
    * \param[in] Wj The second weight \f$ w_j \f$.
    * \param[in] Wk The third weight \f$ w_k \f$.
    * \param[in] Wl The fourth weight \f$ w_l \f$.
    * \param[in] Permutation A reference to a permutation of operators for this part.
    */
    void addMultiterm(ComplexType Coeff, RealType beta,
                      RealType Ei, RealType Ej, RealType Ek, RealType El,
                      RealType Wi, RealType Wj, RealType Wk, RealType Wl);

    /** A difference in energies with magnitude less than this value is treated as zero. */
    static const RealType KroneckerSymbolTolerance = 1e-16;//1e-16;
    /** A difference in energies with magnitude less than this value is treated as zero. */
    static const RealType ReduceResonanceTolerance = 1e-8;//1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account. */
    static const RealType CoefficientTolerance = 1e-16;//1e-16;
    /** A maximum amount of terms put into a list at which the reduceTerms method should be called */
    static const RealType ReduceInvocationThreshold = 1e5;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
    static const RealType MultiTermCoefficientTolerance = 1e-5;//1e-5;

    /** Reduces the number of calculated terms 
    * \param[in] NonResonantTolerance The tolerance for the nonresonant terms cutoff.
    * \param[in] ResonantTolerance    The tolerance for the resonant terms cutoff.
    * \param[in] NonResonantTerms     The list of nonresonant terms.
    * \param[in] ResonantTerms        The list of resonant terms.
    */
    static void reduceTerms(const RealType NonResonantTolerance, const RealType ResonantTolerance, std::list<NonResonantTerm>& NonResonantTerms, std::list<ResonantTerm>& ResonantTerms);

public:
    /** Constructor.
     * \param[in] O1 A reference to a part of the first operator.
     * \param[in] O2 A reference to a part of the second operator.
     * \param[in] O3 A reference to a part of the third operator.
     * \param[in] CX4 A reference to a part of the fourth (creation) operator.
     * \param[in] Hpart1 A reference to the first part of a Hamiltonian.
     * \param[in] Hpart2 A reference to the second part of a Hamiltonian.
     * \param[in] Hpart3 A reference to the third part of a Hamiltonian.
     * \param[in] Hpart4 A reference to the fourth part of a Hamiltonian.
     * \param[in] DMpart1 A reference to the first part of a density matrix.
     * \param[in] DMpart2 A reference to the second part of a density matrix.
     * \param[in] DMpart3 A reference to the third part of a density matrix.
     * \param[in] DMpart4 A reference to the fourth part of a density matrix.
     * \param[in] Permutation A permutation of the operators for this part.
     */
    TwoParticleGFPart(FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation);

    /** Actually computes the part.
     * \param[in] NumberOfMatsubaras: TODO: describes this parameter.
     */
    void compute(long NumberOfMatsubaras);

    /** Returns a contribution to the two-particle Green's function made by this part.
    * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
    * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
    * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
    */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;


    /** Returns the number of resonant terms in the cache. */
    size_t getNumResonantTerms() const;
    /** Returns the number of non-resonant terms in the cache. */
    size_t getNumNonResonantTerms() const;
    
    /** Returns a Permutation3 of the current part */
    const Permutation3& getPermutation();

    /** Return the list of Resonant Terms */
    const std::list<NonResonantTerm>& getNonResonantTerms();
    /** Return the list of NonResonantTerms */
    const std::list<ResonantTerm>& getResonantTerms();
    /** Purges all terms. */
    void clear();
};

struct TwoParticleGFPart::NonResonantTerm{
    /** Coefficient \f$ C \f$. */
    ComplexType Coeff;

    /** Poles \f$ P_1 \f$, \f$ P_2 \f$, \f$ P_3 \f$. */
    RealType Poles[3];
        
    /** Are we using \f$ z_4 \f$ instead of \f$ z_2 \f$ in this term? */
    bool isz4;

    /** Constructor.
    * \param[in] Coeff Numerator of the term.
    * \param[in] P1 Pole P1.
    * \param[in] P2 Pole P2.
    * \param[in] P3 Pole P3.
    * \param[in] isz4 Are we using \f$ z_4 \f$ instead of \f$ z_2 \f$ in this term?
    */
    NonResonantTerm(ComplexType Coeff, RealType P1, RealType P2, RealType P3, bool isz4);

    /** Returns a contribution to the two-particle Green's function made by this term.
    * \param[in] z1 Complex frequency \f$ z_1 \f$.
    * \param[in] z2 Complex frequency \f$ z_2 \f$.
    * \param[in] z3 Complex frequency \f$ z_3 \f$.
    */
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;

    /** This operator add a non-resonant term to this one.
    * It does not check the similarity of the terms! 
    * \param[in] AnotherTerm Another term to add to this.
    */
    NonResonantTerm& operator+=(const NonResonantTerm& AnotherTerm);
        
    /** Returns true if another term is similar to this
     * (sum of the terms is again a correct non-resonant term).
    */
    bool isSimilarTo(const NonResonantTerm& AnotherTerm) const;
    };
  
struct TwoParticleGFPart::ResonantTerm {

    /** Coefficient \f$ R \f$. */
    ComplexType ResCoeff;
    /** Coefficient \f$ N \f$. */
    ComplexType NonResCoeff;

    /** Poles \f$ P_1 \f$, \f$ P_2 \f$, \f$ P_3 \f$. */
    RealType Poles[3];
        
    /** Are we using \f$ \delta(z_1+z_2-P_1-P_2) \f$ resonance condition?
     Otherwise we are using \f$ \delta(z_2+z_3-P_2-P_3) \f$. */
    bool isz1z2;

    /** Constructor.
     * \param[in] ResCoeff Numerator of the term for a resonant case.
     * \param[in] NonResCoeff Numerator of the term for a non-resonant case.
     * \param[in] P1 Pole P1.
     * \param[in] P2 Pole P2.
     * \param[in] P3 Pole P3.
     * \param[in] isz1z2 Are we using \f$ \delta(z_1+z_2-P_1-P_2) \f$ resonance condition?
     */
    ResonantTerm(ComplexType ResCoeff, ComplexType NonResCoeff, RealType P1, RealType P2, RealType P3, bool isz1z2);

    /** Returns a contribution to the two-particle Green's function made by this term.
    * \param[in] z1 Complex frequency \f$ z_1 \f$.
    * \param[in] z2 Complex frequency \f$ z_2 \f$.
    * \param[in] z3 Complex frequency \f$ z_3 \f$.
    */
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;
        
    /** This operator add a non-resonant term to this one.
    * It does not check the similarity of the terms! 
    * \param[in] AnotherTerm Another term to add to this.
    */
    ResonantTerm& operator+=(const ResonantTerm& AnotherTerm);
        
    /** Returns true if another term is similar to this
     * (sum of the terms is again a correct non-resonant term).
    */
    bool isSimilarTo(const ResonantTerm& AnotherTerm) const;
};

inline
ComplexType TwoParticleGFPart::NonResonantTerm::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const
{
    return isz4 ?   Coeff / ((z1-Poles[0])*(z1+z2+z3-Poles[0]-Poles[1]-Poles[2])*(z3-Poles[2])) :
                    Coeff / ((z1-Poles[0])*(z2-Poles[1])*(z3-Poles[2]));
}

inline
ComplexType TwoParticleGFPart::ResonantTerm::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const
{
    ComplexType Diff;
    if(isz1z2){
        Diff = z1 + z2 - Poles[0] - Poles[1];
        return (abs(Diff) < KroneckerSymbolTolerance ? ResCoeff : (NonResCoeff/Diff) )
                /((z1-Poles[0])*(z3-Poles[2]));
    } else {
        Diff = z2 + z3 - Poles[1] - Poles[2];
        return (abs(Diff) < KroneckerSymbolTolerance ? ResCoeff : (NonResCoeff/Diff) )
                /((z1-Poles[0])*(z3-Poles[2]));
    }
}

#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGFPART_H
