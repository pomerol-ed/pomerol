//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/HamiltonianPart.hpp
/// \brief Storage and diagonalization of a Hamiltonian matrix block.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_HPP
#define POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_HPP

#include "ComputableObject.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <iostream>
#include <memory>
#include <type_traits>

namespace Pomerol {

/// \addtogroup ED
///@{

/// \brief Part of a Hamiltonian of a quantum system.
///
/// This class stores and diagonalizes a single block of the Hamiltonian matrix, which corresponds to a single
/// invariant subspace of the Hamiltonian.
class HamiltonianPart : public ComputableObject {

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;

    /// Index of the block (invariant subspace) this part corresponds to.
    BlockNumber Block;

    /// Whether the stored matrix is complex-valued.
    bool Complex;

    /// A type-erased pointer to the respective real/complex-valued \p libcommute::loperator object.
    void const* HOp;

    /// The type-erased real/complex matrix of this block of the Hamiltonian.
    std::shared_ptr<void> HMatrix = nullptr;

    /// Eigenvalues of this block.
    RealVectorType Eigenvalues;

    friend class Hamiltonian;

public:
    /// Constructor.
    /// \tparam ScalarType Scalar type (either double or std::complex<double>) of the linear operator \p HOp.
    /// \param[in] HOp The linear operator object corresponding to the Hamiltonian.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] Block Index of the block (invariant subspace) this part corresponds to.
    template <typename ScalarType>
    HamiltonianPart(LOperatorType<ScalarType> const& HOp, StatesClassification const& S, BlockNumber Block)
        : S(S), Block(Block), Complex(std::is_same<ScalarType, ComplexType>::value), HOp(&HOp) {}

    /// Fill the matrix with elements.
    void prepare();

    /// Diagonalize the matrix.
    /// \pre \ref prepare() has been called.
    void compute();

    /// Discard all eigenvalues exceeding a given cutoff and truncate the size of the diagonalized
    /// matrix accordingly.
    /// \param[in] Cutoff Maximum allowed value of the energy.
    /// \pre \ref compute() has been called.
    bool reduce(RealType Cutoff);

    /// Is this object storing a complex-valued matrix?
    bool isComplex() const { return Complex; }

    /// Return the index of the block (invariant subspace) this part corresponds to.
    BlockNumber getBlockNumber() const { return Block; }

    /// Return dimension of the respective invariant subspace.
    InnerQuantumState getSize() const;

    /// Access eigenvalues of the matrix.
    /// \pre \ref compute() has been called.
    RealVectorType const& getEigenValues() const;

    /// Access a single eigenvalue.
    /// \param[in] State Index of the eigenvalue.
    /// \pre \ref compute() has been called.
    RealType getEigenValue(InnerQuantumState State) const;

    /// Return a constant reference to the stored matrix.
    /// \tparam Complex Request a reference to a complex-valued matrix.
    /// \pre \ref prepare() has been called.
    /// \pre The compile-time value of \p Complex must agree with the result of \ref isComplex().
    template <bool Complex> MatrixType<Complex> const& getMatrix() const;
    /// Return a reference to the stored matrix.
    /// \tparam Complex Request a reference to a complex-valued matrix.
    /// \pre \ref prepare() has been called.
    /// \pre The compile-time value of \p Complex must agree with the result of \ref isComplex().
    template <bool Complex> MatrixType<Complex>& getMatrix();

    /// Return the lowest eigenvalue.
    /// \pre \ref compute() has been called.
    RealType getMinimumEigenvalue() const;

    /// Return a single eigenstate.
    /// \tparam Complex Request a reference to a complex-valued eigenvector.
    /// \param[in] State Index of the eigenstate.
    /// \pre The compile-time value of \p Complex must agree with the result of \ref isComplex().
    template <bool Complex> VectorType<Complex> getEigenState(InnerQuantumState State) const;

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] part HamiltonianPart to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, HamiltonianPart const& part) {
        if(part.isComplex())
            os << part.getMatrix<true>() << std::endl;
        else
            os << part.getMatrix<false>() << std::endl;

        return os;
    }

private:
    // Implementation details
    template <bool C> void initHMatrix();
    template <bool C> void prepareImpl();
    template <bool C> void computeImpl();

    void checkComputed() const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_HPP
