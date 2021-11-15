//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/Hamiltonian.h
** \brief A storage of the matrix elements of the hamiltonian in Fock basis, provides eigenvalues and eigenfunctions
**
** \author Andrey Antipov(Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIAN_HPP
#define POMEROL_INCLUDE_POMEROL_HAMILTONIAN_HPP

#include "ComputableObject.hpp"
#include "HamiltonianPart.hpp"
#include "HilbertSpace.hpp"
#include "Misc.hpp"
#include "Operators.hpp"
#include "StatesClassification.hpp"

#include "mpi_dispatcher/misc.hpp"

#include <cmath>
#include <type_traits>
#include <vector>

namespace Pomerol {

/** This class represents a Hamiltonian, written as a matrix of matrix elements in a Fock basis.
 * It is a container for several hamiltonian parts, each for single defined QuantumNumbers and a corresponding BlockNumber.
 * It provides eigenvalues and eigenfunctions of any of its parts once they are obtained within its parts.
 * The diagonalization and entering routines are done inside of HamiltonianPart instances.
 */
class Hamiltonian : public ComputableObject {
    bool Complex = {};

    /** Array of pointers to the Hamiltonian Parts */
    std::vector<HamiltonianPart> parts;
    /** A reference to the StatesClassification object. */
    StatesClassification const& S;
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy = -HUGE_VAL;

public:
    /** Constructor. */
    explicit Hamiltonian(StatesClassification const& S) : S(S) {}

    template <typename ScalarType, typename... IndexTypes>
    void prepare(Operators::expression<ScalarType, IndexTypes...> const& H,
                 HilbertSpace<IndexTypes...> const& HS,
                 MPI_Comm const& comm = MPI_COMM_WORLD);
    void compute(MPI_Comm const& comm = MPI_COMM_WORLD);
    void reduce(RealType Cutoff);

    bool isComplex() const { return Complex; }

    HamiltonianPart const& getPart(BlockNumber Block) const { return parts[Block]; }

    InnerQuantumState getBlockSize(BlockNumber Block) const;

    RealType getEigenValue(unsigned long state) const;
    RealVectorType const& getEigenValues(BlockNumber Block) const;
    RealVectorType getEigenValues() const;
    RealType getGroundEnergy() const { return GroundEnergy; }

private:
    void computeGroundEnergy();

    template <bool C> void prepareImpl(LOperatorTypeRC<C> const& HOp, const MPI_Comm& comm);
    template <bool C> void computeImpl(MPI_Comm const& comm);
};

template <typename ScalarType, typename... IndexTypes>
void Hamiltonian::prepare(Operators::expression<ScalarType, IndexTypes...> const& H,
                          HilbertSpace<IndexTypes...> const& HS,
                          MPI_Comm const& comm) {

    if(getStatus() >= Prepared)
        return;

    Complex = std::is_same<ScalarType, ComplexType>::value;
    LOperatorType<ScalarType> HOp(H, HS.getFullHilbertSpace());
    prepareImpl<std::is_same<ScalarType, ComplexType>::value>(HOp, comm);

    setStatus(Prepared);
}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIAN_HPP
