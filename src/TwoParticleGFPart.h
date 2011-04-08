/** \file src/TwoParticleGFPart.h
** \brief Part of a two-particle Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef ____DEFINE_2PGF_PART____
#define ____DEFINE_2PGF_PART____

#include <list>

#include "config.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "FieldOperator.h"
#include "DensityMatrixPart.h"

/** This class represents a part of a two-particle Green's function.
 * Every part describes one 'world stripe' of four operators.
 */
class TwoParticleGFPart {

public:

    /** A summation over all matrix elements may be performed (currently) in two different orders.
     * So there are two implemented methods of computation.
     */
    enum ComputationMethod {ChasingIndices1, ChasingIndices2};

    /** Every term is a sum of fractions:
     * \f[
     * \frac{1}{(z_1-P_1)(z_3-P_3)}
     *         \left(\frac{C_4}{z_1+z_2+z_3-P_1-P_2-P_3} - \frac{C_2}{z_2-P_2} \right. +
     * \f]
     * \f[     \left.
     *         + R_{12}\delta(z_1+z_2-P_1-P_2)
     *         + N_{12}\frac{1 - \delta(z_1+z_2-P_1-P_2)}{z_1+z_2-P_1-P_2}
     *         + R_{23}\delta(z_2+z_3-P_2-P_3)
     *         + N_{23}\frac{1 - \delta(z_2+z_3-P_2-P_3)}{z_2+z_3-P_2-P_3}
     *         \right)
     * \f]
     *
     * In fact this is a slightly rewritten form of an equation for \f$ \phi \f$ from
     * <em>H. Hafermann et al 2009 EPL 85 27007</em>.
     */
    struct TwoParticleGFTerm{
        /** The index of the frequency which will substitute \f$ z_1 \f$. */
        size_t z1;
        /** The index of the frequency which will substitute \f$ z_2 \f$. */
        size_t z2;
        /** The index of the frequency which will substitute \f$ z_3 \f$. */
        size_t z3;

        /** Coefficient \f$ C_2.\f$ */
        ComplexType CoeffZ2;
        /** Coefficient \f$ C_4.\f$ */
        ComplexType CoeffZ4;
        /** Coefficient \f$ R_{12}.\f$ */
        ComplexType CoeffZ1Z2Res;
        /** Coefficient \f$ N_{12}.\f$ */
        ComplexType CoeffZ1Z2NonRes;
        /** Coefficient \f$ R_{23}.\f$ */
        ComplexType CoeffZ2Z3Res;
        /** Coefficient \f$ N_{23}.\f$ */
        ComplexType CoeffZ2Z3NonRes;

        /** Array of poles \f$ P_1 \f$, \f$ P_2 \f$, \f$ P_3 \f$. */
        ComplexType Poles[3];

        /** Constructor.
         *  It calculates all the member coefficients and poles as follows:
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
        TwoParticleGFTerm(ComplexType Coeff, RealType beta,
                         RealType Ei, RealType Ej, RealType Ek, RealType El,
                         RealType Wi, RealType Wj, RealType Wk, RealType Wl,
                         Permutation3& Permutation);

        /** Returns a contribution to the two-particle Green's function made by this term.
        * \param[in] Frequency1 Complex frequency \f$ i\omega_1 \f$ to substitute into this term.
        * \param[in] Frequency2 Complex frequency \f$ i\omega_2 \f$ to substitute into this term.
        * \param[in] Frequency3 Complex frequency \f$ i\omega_3 \f$ to substitute into this term.
        */
        ComplexType operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const;

        /** A matrix element with magnitude less than this value is treated as zero. */
        static const RealType MatrixElementTolerance = 1e-3;//1e-8;//TERM_MATRIX_ELEMENT_TOLERANCE;
        /** A difference in energies with magnitude less than this value is treated as zero. */
        static const RealType ResonanceTolerance = 1e-16;//TERM_RESONANCE_TOLERANCE ;
        /** Should we take a product of matrix elements into account?
         * \param[in] MatrixElementProd A product of four matrix elements.
         */
        static bool IsRelevantMatrixElement(const ComplexType &MatrixElementProd);
        bool hasRelevantNumerator(const RealType Tolerance);
    };

    /**
     * A miniclass to store value of Chi over Matsubara frequencies. Stores data in a (OMEGA,nu,nu'), where OMEGA=w1+w2 - bosonic frequency
     * and nu=w1, nu'=w4
     *
     * TODO: document this class in depth.
     */
    class MatsubaraContainer{

        ComplexType MatsubaraSpacing;
        long NumberOfMatsubaras;
        std::vector<MatrixType> Data;
        std::vector<long> FermionicFirstIndex;
    public:
        MatsubaraContainer(RealType beta);
        void prepare(long NumberOfMatsubaras);
        ComplexType& operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
        inline void fill(std::list<TwoParticleGFTerm> &Terms);
        MatsubaraContainer& operator+=(const MatsubaraContainer& rhs);
        void clear();
    };

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

    /** A list of pointers to parts. */
    std::list<TwoParticleGFTerm> Terms;

    /** TODO: document this */
    MatsubaraContainer *Storage;

    /** Computes the part using the ChasingIndices1 method.
     * \param[in] NumberOfMatsubaras: TODO: describe this parameter.
     */
    void computeChasing1(long NumberOfMatsubaras);
    /** Computes the part using the ChasingIndices2 method.
     * \param[in] NumberOfMatsubaras: TODO: describe this parameter.
     */
    void computeChasing2(long NumberOfMatsubaras);

    /** Reduces the number of calculated terms */
    void reduceTerms();

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
     * \param[in] method Method of computation to use.
     */
    void compute(long NumberOfMatsubaras, ComputationMethod method = ChasingIndices2);
    /** Purges all terms. */
    void clear();
     /** Returns a contribution to the two-particle Green's function made by this part.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    /** TODO: describe this method. */
    const MatsubaraContainer& getMatsubaraContainer();
};

#endif // endif :: #ifndef ____DEFINE_2PGF_PART____
