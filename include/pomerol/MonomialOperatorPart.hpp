//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/MonomialOperatorPart.hpp
/// \brief Storage for a matrix block of an operator that is a product of creation/annihilation operators.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_MONOMIALOPERATORPART_HPP
#define POMEROL_INCLUDE_MONOMIALOPERATORPART_HPP

#include "HamiltonianPart.hpp"
#include "HilbertSpace.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <memory>
#include <ostream>
#include <type_traits>

namespace Pomerol {

/// \addtogroup ED
///@{

/// \brief Part of a monomial quantum operator.
///
/// This class stores a matrix block of an operator \f$\hat M\f$, which is a monomial, i.e. a product
/// of fermionic and/or bosonic creation/annihilation operators. The matrix is computed in the eigenbasis
/// of the Hamiltonian \f$\hat H\f$ and connects one of its invariant subspaces (right subspace)
/// to another one (left subspace),
/// \f[
///   \langle {\rm left}|\hat M|{\rm right}\rangle.
/// \f]
class MonomialOperatorPart : public ComputableObject {
    friend class FieldOperatorContainer;

private:
    /// Whether the following \p libcommute::loperator object is complex-valued.
    bool MOpComplex;
    /// A type-erased pointer to the respective real/complex-valued \p libcommute::loperator object.
    void const* MOp;

    /// Whether the stored matrices are complex-valued.
    bool Complex;

    /// Information about invariant subspaces of the Hamiltonian.
    StatesClassification const& S;
    /// Diagonal block of the Hamiltonian corresponding to the right invariant subspace.
    HamiltonianPart const& HFrom;
    /// Diagonal block of the Hamiltonian corresponding to the left invariant subspace.
    HamiltonianPart const& HTo;

protected:
    /// Type-erased real/complex sparse matrix \f$\langle {\rm left}|\hat M|{\rm right}\rangle\f$
    /// stored in the row-major order.
    std::shared_ptr<void> elementsRowMajor;
    /// Type-erased real/complex sparse matrix \f$\langle {\rm left}|\hat M|{\rm right}\rangle\f$
    /// stored in the column-major order.
    std::shared_ptr<void> elementsColMajor;

public:
    /// Constructor.
    /// \tparam ScalarType Scalar type (either double or std::complex<double>) of the linear operator \p MOp.
    /// \param[in] MOp The linear operator object corresponding to the monomial operator \f$\hat M\f$.
    /// \param[in] S Information about invariant subspaces of the Hamiltonian.
    /// \param[in] HFrom Diagonal block of the Hamiltonian corresponding to the right invariant subspace.
    /// \param[in] HTo Diagonal block of the Hamiltonian corresponding to the left invariant subspace.
    template <typename ScalarType>
    MonomialOperatorPart(LOperatorType<ScalarType> const& MOp,
                         StatesClassification const& S,
                         HamiltonianPart const& HFrom,
                         HamiltonianPart const& HTo)
        : MOpComplex(std::is_same<ScalarType, ComplexType>::value),
          MOp(&MOp),
          Complex(MOpComplex || HFrom.isComplex() || HTo.isComplex()),
          S(S),
          HFrom(HFrom),
          HTo(HTo) {}

    /// Compute and store all matrix elements of \f$\hat M\f$ in the eigenbasis of the Hamiltonian.
    /// \param[in] Tolerance Matrix elements with the absolute value equal or below this threshold
    ///                      are considered negligible.
    void compute(RealType Tolerance);

    /// Reset the stored sparse matrices to those obtained from
    /// \f$\langle {\rm right}|\hat M^\dagger|{\rm left}\rangle\f$.
    /// \param[in] part Monomial operator part \f$\langle {\rm right}|\hat M^\dagger|{\rm left}\rangle\f$.
    void setFromAdjoint(MonomialOperatorPart const& part);

    /// Is this object storing a complex-valued sparse matrices?
    bool isComplex() const { return Complex; }

    /// Return a reference to the stored row-major sparse matrix.
    /// \tparam C Request a reference to the complex-valued matrix.
    /// \pre The compile-time value of \p C must agree with the result of \ref isComplex().
    template <bool C> RowMajorMatrixType<C>& getRowMajorValue();
    /// Return a constant reference to the stored row-major sparse matrix.
    /// \tparam C Request a reference to the complex-valued matrix.
    /// \pre The compile-time value of \p C must agree with the result of \ref isComplex().
    template <bool C> RowMajorMatrixType<C> const& getRowMajorValue() const;
    /// Return a reference to the stored column-major sparse matrix.
    /// \tparam C Request a reference to the complex-valued matrix.
    /// \pre The compile-time value of \p C must agree with the result of \ref isComplex().
    template <bool C> ColMajorMatrixType<C>& getColMajorValue();
    /// Return a constant reference to the stored column-major sparse matrix.
    /// \tparam C Request a reference to the complex-valued matrix.
    /// \pre The compile-time value of \p C must agree with the result of \ref isComplex().
    template <bool C> ColMajorMatrixType<C> const& getColMajorValue() const;

    /// Return the index of the right invariant subspace.
    BlockNumber getRightIndex() const { return HFrom.getBlockNumber(); }
    /// Return the index of the left invariant subspace.
    BlockNumber getLeftIndex() const { return HTo.getBlockNumber(); }

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] part \ref MonomialOperatorPart to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, MonomialOperatorPart const& part) {
        if(part.isComplex())
            part.streamOutputImpl<true>(os);
        else
            part.streamOutputImpl<false>(os);
        return os;
    }

private:
    // Implementation details
    template <bool C, bool HC> void computeImpl(RealType Tolerance);
    template <bool C> void streamOutputImpl(std::ostream& os) const;
};

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_MONOMIALOPERATORPART_HPP
