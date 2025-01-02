//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Hamiltonian.hpp
/// \brief Storage and diagonalization of a Hamiltonian matrix.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

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

/// \addtogroup ED
///@{

/// \brief Hamiltonian of a quantum system.
///
/// This class represents a Hamiltonian as a block-diagonal matrix with blocks corresponding
/// to distinct invariant subspaces. The blocks are stored as a list of
/// \ref HamiltonianPart objects. The main purpose of this class is MPI-parallelized
/// diagonalization of the entire Hamiltonian matrix.
class Hamiltonian : public ComputableObject {

    /// Whether the Hamiltonian is complex-valued.
    bool Complex = {};

    /// List of parts (diagonal matrix blocks).
    std::vector<HamiltonianPart> parts;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;

    /// The ground state energy.
    RealType GroundEnergy = -HUGE_VAL;

public:
    /// Constructor.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    explicit Hamiltonian(StatesClassification const& S) : S(S) {}

    /// Fill matrices of all diagonal blocks in parallel.
    /// \tparam ScalarType Scalar type (either double or std::complex<double>) of the expression \p H.
    /// \tparam IndexTypes Types of indices carried by operators in the expression \p H.
    /// \param[in] H Expression of the Hamiltonian.
    /// \param[in] HS Hilbert space.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    template <typename ScalarType, typename... IndexTypes>
    void prepare(Operators::expression<ScalarType, IndexTypes...> const& H,
                 HilbertSpace<IndexTypes...> const& HS,
                 MPI_Comm const& comm = MPI_COMM_WORLD);

    /// Diagonalize matrices of all diagonal blocks in parallel.
    /// \param[in] comm MPI communicator used to parallelize the computation.
    /// \pre \ref prepare() has been called.
    void compute(MPI_Comm const& comm = MPI_COMM_WORLD);

    /// Discard all eigenvalues exceeding a given cutoff and truncate the size of all diagonalized
    /// blocks accordingly.
    /// \param[in] Cutoff Maximum allowed excitation energy (energy level calculated w.r.t. the ground state energy).
    /// \pre \ref compute() has been called.
    void reduce(RealType Cutoff);

    /// Is the Hamiltonian a complex-valued matrix?
    bool isComplex() const { return Complex; }

    /// Access a part (diagonal block) of the Hamiltonian.
    /// \param[in] Block Index of the diagonal block.
    /// \pre \ref prepare() has been called.
    HamiltonianPart const& getPart(BlockNumber Block) const { return parts[Block]; }

    /// Return size of a part (dimension of a diagonal block).
    /// \param[in] Block Index of the diagonal block.
    /// \pre \ref prepare() has been called.
    InnerQuantumState getBlockSize(BlockNumber Block) const;

    /// Return a single eigenvalue of the Hamiltonian.
    /// \param[in] state Index of the eigenvalue within the full diagonalized matrix of the Hamiltonian.
    /// \pre \ref compute() has been called.
    RealType getEigenValue(QuantumState state) const;

    /// Return a list of eigenvalues of the Hamiltonian within a block.
    /// \param[in] Block Index of the diagonal block.
    /// \pre \ref compute() has been called.
    RealVectorType const& getEigenValues(BlockNumber Block) const;

    /// Return a list of all eigenvalues of the Hamiltonian.
    /// \pre \ref compute() has been called.
    RealVectorType getEigenValues() const;

    /// Return the ground state energy.
    /// \pre \ref compute() has been called.
    RealType getGroundEnergy() const { return GroundEnergy; }

private:
    // Implementation details
    void computeGroundEnergy();

    // cppcheck-suppress unusedPrivateFunction
    template <bool C> void prepareImpl(LOperatorTypeRC<C> const& HOp, MPI_Comm const& comm);
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

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIAN_HPP
