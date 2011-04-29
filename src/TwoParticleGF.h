/** \file src/TwoParticleGF.h
** \brief Two-particle Green's function in the Matsubara representation.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#ifndef __INCLUDE_TWOPARTICLEGF_H
#define __INCLUDE_TWOPARTICLEGF_H

#include <sstream>
#include <fstream>
#include <iomanip>

#include"Misc.h"
#include"ComputableObject.h"
#include"FourIndexObject.h"
#include"output.h"
#include"StatesClassification.h"
#include"FieldOperator.h"
#include"DensityMatrix.h"
#include"TwoParticleGFPart.h"

/** This class represents a thermal Green's function in the Matsubara representation.
 *
 * Exact definition:
 *
 * \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
 *   -\int_0^\beta
 *      \langle\mathbf{T} c_i(\tau_1)c_j(\tau_2)c^+_k(\tau_3)c^+_l(0) \rangle
 *      \exp(i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3)
 *    d\tau_1 d\tau_2 d\tau_3
 * \f].
 *
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). Every part corresponds to a 'world-stripe', a sequence of 4
 * matrix blocks.
 */
class TwoParticleGF : public ComputableObject, public FourIndexSingleObject, public Thermal {

    /** A reference to a states classification object. */
    StatesClassification& S;
    /** A reference to a Hamiltonian. */
    Hamiltonian& H;
    /** A reference to the first annihilation operator. */
    AnnihilationOperator& C1;
    /** A reference to the second annihilation operator. */
    AnnihilationOperator& C2;
    /** A reference to the first creation operator. */
    CreationOperator& CX3;
    /** A reference to the second creation operator. */
    CreationOperator& CX4;
    /** A reference to a density matrix. */
    DensityMatrix& DM;

    /** A list of pointers to parts. */
    std::list<TwoParticleGFPart*> parts;

    std::list<TwoParticleGFPart::ResonantTerm> ResonantTerms[6];
    std::list<TwoParticleGFPart::NonResonantTerm> NonResonantTerms[6];

    /** A flag to determine whether this GF is identical to zero */
    bool vanish;

    /** Extracts a part of the operator standing at a specified position in a given permutation.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] LeftIndex A left block index referring to the part needed.
     */
    FieldOperatorPart& OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex);
    /** Chooses an operator standing at a specified position in a given permutation and
     * returns a left block index corresponding to the right block index. May return ERROR_BLOCK_NUMBER if
     * the operator does not have such a (non-zero) block.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] RightIndex A right block index.
     */
    BlockNumber getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex);
    /** Chooses an operator standing at a specified position in a given permutation and
     * returns a right block index corresponding to the left block index. May return ERROR_BLOCK_NUMBER if
     * the operator does not have such a (non-zero) block.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] LeftIndex A left block index.
     */
    BlockNumber getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex); //!< return right index of an operator at current position for a current permutation

    /** TODO: need a description. */
    MatsubaraContainer *Storage;

public:
    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] C1 A reference to the first annihilation operator.
     * \param[in] C2 A reference to the second annihilation operator.
     * \param[in] CX3 A reference to the first creation operator.
     * \param[in] CX4 A reference to the second creation operator.
     * \param[in] DM A reference to a density matrix.
     */
    TwoParticleGF(StatesClassification& S, Hamiltonian& H,
            AnnihilationOperator& C1, AnnihilationOperator& C2,
            CreationOperator& CX3, CreationOperator& CX4,
            DensityMatrix& DM);
    /** Destructor. */
    ~TwoParticleGF();

    /** Chooses relevant parts of C1, C2, CX3 and CX4 and allocates resources for the parts. */
    void prepare(void);
    /** Actually computes the parts.
     * \param[in] NumberOfMatsubaras: TODO: need a description.
     */
    void compute(long NumberOfMatsubaras);

    /** Returns the 'bit' (index) of one of operators C1, C2, CX3 or CX4.
     * \param[in] Position Zero-based number of the operator to use.
     */
    ParticleIndex getIndex(size_t Position) const;

    /** Returns the value of the two-particle Green's function calculated at given frequencies.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    /** Returns true, if GF is identical to zero */
    bool vanishes();

    /** Returns the total number of resonant terms for all parts. */
    size_t getNumResonantTerms() const;
    /** Returns the totla number of non-resonant terms for all parts. */
    size_t getNumNonResonantTerms() const;
    
    /** Returns the number of current permutation in permutations3 */
    unsigned short getPermutationNumber(const Permutation3& in);
};

#endif // endif :: #ifndef __INCLUDE_TWOPARTICLEGF_H
