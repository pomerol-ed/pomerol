//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/TwoParticleGF.h
** \brief Two-particle Green's function in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_TWOPARTICLEGF_H
#define POMEROL_INCLUDE_TWOPARTICLEGF_H

#include "ComputableObject.hpp"
#include "DensityMatrix.hpp"
#include "Hamiltonian.hpp"
#include "Misc.hpp"
#include "MonomialOperator.hpp"
#include "StatesClassification.hpp"
#include "Thermal.hpp"
#include "TwoParticleGFPart.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <cstddef>
#include <numeric>
#include <tuple>
#include <vector>

namespace Pomerol {

using FreqTuple = std::tuple<ComplexType, ComplexType, ComplexType>;
using FreqVec = std::vector<FreqTuple>;

/** This class represents a thermal Green's function in the Matsubara representation.
 *
 * Exact definition:
 *
 * \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
 *   \int_0^\beta
 *     \langle\mathbf{T} c_i(\tau_1)c_j(\tau_2)c^+_k(\tau_3)c^+_l(0) \rangle
 *     \exp(i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3)
 *   d\tau_1 d\tau_2 d\tau_3
 * \f].
 *
 * It is actually a container class for a collection of parts (most of real calculations
 * take place inside the parts). Every part corresponds to a 'world-stripe', a sequence of 4
 * matrix blocks.
 */
class TwoParticleGF : public Thermal, public ComputableObject {

    friend class TwoParticleGFContainer;

    /** A reference to a states classification object. */
    StatesClassification const& S;
    /** A reference to a Hamiltonian. */
    Hamiltonian const& H;
    /** A reference to the first annihilation operator. */
    AnnihilationOperator const& C1;
    /** A reference to the second annihilation operator. */
    AnnihilationOperator const& C2;
    /** A reference to the first creation operator. */
    CreationOperator const& CX3;
    /** A reference to the second creation operator. */
    CreationOperator const& CX4;
    /** A reference to a density matrix. */
    DensityMatrix const& DM;

    /** A list of pointers to parts. */
    std::vector<TwoParticleGFPart> parts;

protected:
    /** A flag to determine whether this GF is identical to zero */
    bool Vanishing = true;

    /** Extracts a part of the operator standing at a specified position in a given permutation.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] LeftIndex A left block index referring to the part needed.
     */
    MonomialOperatorPart const&
    OperatorPartAtPosition(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber LeftIndex) const;
    /** Chooses an operator standing at a specified position in a given permutation and
     * returns a left block index corresponding to the right block index. May return INVALID_BLOCK_NUMBER if
     * the operator does not have such a (non-zero) block.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] RightIndex A right block index.
     */
    BlockNumber getLeftIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber RightIndex) const;
    /** Chooses an operator standing at a specified position in a given permutation and
     * returns a right block index corresponding to the left block index. May return INVALID_BLOCK_NUMBER if
     * the operator does not have such a (non-zero) block.
     * \param[in] PermutationNumber The number of the permutation.
     * \param[in] OperatorPosition The number of the position of the operator.
     * \param[in] LeftIndex A left block index.
     */
    BlockNumber getRightIndex(std::size_t PermutationNumber, std::size_t OperatorPosition, BlockNumber LeftIndex)
        const; //!< return right index of an operator at current position for a current permutation

public:
    /** A difference in energies with magnitude less than this value is treated as zero. default = 1e-8. */
    RealType ReduceResonanceTolerance = 1e-8;
    /** Minimal magnitude of the coefficient of a term to take it into account. default = 1e-16. */
    RealType CoefficientTolerance = 1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. default = 1e-5. */
    RealType MultiTermCoefficientTolerance = 1e-5;

    /** Constructor.
     * \param[in] S A reference to a states classification object.
     * \param[in] H A reference to a Hamiltonian.
     * \param[in] C1 A reference to the first annihilation operator.
     * \param[in] C2 A reference to the second annihilation operator.
     * \param[in] CX3 A reference to the first creation operator.
     * \param[in] CX4 A reference to the second creation operator.
     * \param[in] DM A reference to a density matrix.
     */
    TwoParticleGF(StatesClassification const& S,
                  Hamiltonian const& H,
                  AnnihilationOperator const& C1,
                  AnnihilationOperator const& C2,
                  CreationOperator const& CX3,
                  CreationOperator const& CX4,
                  DensityMatrix const& DM);

    /** Chooses relevant parts of C1, C2, CX3 and CX4 and allocates resources for the parts. */
    void prepare();

    /** Actually computes the parts and fill the internal cache of precomputed values.
     * \param[in] NumberOfMatsubaras Number of positive Matsubara frequencies.
     */
    std::vector<ComplexType>
    compute(bool clear = false, FreqVec const& freqs = {}, MPI_Comm const& comm = MPI_COMM_WORLD);

    /** Returns the 'bit' (index) of one of operators C1, C2, CX3 or CX4.
     * \param[in] Position Zero-based number of the operator to use.
     */
    ParticleIndex getIndex(std::size_t Position) const;

    /** Returns the value of the Green's function calculated at a given frequency (ignores precomputed values).
    * \param[in] z1 Frequency 1
    * \param[in] z2 Frequency 2
    * \param[in] z3 Frequency 3
    */
    ComplexType operator()(ComplexType z1, ComplexType z2, ComplexType z3) const;
    /** Returns the value of the two-particle Green's function calculated at given frequencies.
     * \param[in] MatsubaraNumber1 Number of the first Matsubara frequency.
     * \param[in] MatsubaraNumber2 Number of the second Matsubara frequency.
     * \param[in] MatsubaraNumber3 Number of the third Matsubara frequency.
     */
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;

    /** Returns true, if GF is identical to zero */
    bool isVanishing() const { return Vanishing; }

    /** Returns the number of current permutation in permutations3 */
    unsigned short getPermutationNumber(Permutation3 const& in);
};

inline ComplexType TwoParticleGF::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const {
    if(Vanishing)
        return 0;
    else {
        return std::accumulate(parts.begin(),
                               parts.end(),
                               ComplexType(0),
                               [z1, z2, z3](ComplexType s, TwoParticleGFPart const& p) { return s + p(z1, z2, z3); });
    }
}

inline ComplexType
TwoParticleGF::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const {
    return (*this)(MatsubaraSpacing * RealType(2 * MatsubaraNumber1 + 1),
                   MatsubaraSpacing * RealType(2 * MatsubaraNumber2 + 1),
                   MatsubaraSpacing * RealType(2 * MatsubaraNumber3 + 1));
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGF_H
