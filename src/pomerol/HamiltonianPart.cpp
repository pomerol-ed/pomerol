//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/HamiltonianPart.cpp
/// \brief Storage and diagonalization of a Hamiltonian matrix block (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include "pomerol/HamiltonianPart.hpp"

// clang-format off
#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/loperator/mapped_basis_view.hpp>
// clang-format on

#include <Eigen/Eigenvalues>

#include <cassert>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace Pomerol {

//
// class HamiltonianPart
//

template <bool C> void HamiltonianPart::initHMatrix() {
    InnerQuantumState BlockSize = S.getBlockSize(Block);
    HMatrix = std::make_shared<MatrixType<C>>(BlockSize, BlockSize);
}
template void HamiltonianPart::initHMatrix<true>();
template void HamiltonianPart::initHMatrix<false>();

void HamiltonianPart::prepare() {
    if(getStatus() >= Prepared)
        return;

    if(isComplex())
        prepareImpl<true>();
    else
        prepareImpl<false>();

    setStatus(Prepared);
}

template <bool C> void HamiltonianPart::prepareImpl() {
    initHMatrix<C>();

    auto const& HOp_ = *static_cast<LOperatorTypeRC<C> const*>(HOp);
    auto& HMatrix_ = getMatrix<C>();

    auto mapper = libcommute::basis_mapper(S.getFockStates(Block));

    auto BlockSize = S.getBlockSize(Block);
    VectorType<C> ket = VectorType<C>::Zero(BlockSize);
    auto ket_view = mapper.make_const_view(ket);

    for(InnerQuantumState st = 0; st < BlockSize; ++st) {
        auto bra_view = mapper.make_view(HMatrix_.col(st));
        ket(st) = 1.0;
        HOp_(ket_view, bra_view);
        ket(st) = .0; // cppcheck-suppress redundantAssignment
    }

    assert((HMatrix_.adjoint() - HMatrix_).array().abs().maxCoeff() < 100 * std::numeric_limits<RealType>::epsilon());
}

void HamiltonianPart::compute() {
    if(getStatus() >= Computed)
        return;

    if(isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();

    setStatus(Computed);
}

template <bool C> void HamiltonianPart::computeImpl() {
    auto& HMatrix_ = getMatrix<C>();
    if(HMatrix_.rows() == 1) {
        assert(std::abs(HMatrix_(0, 0) - std::real(HMatrix_(0, 0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << std::real(HMatrix_(0, 0));
        HMatrix_(0, 0) = 1;
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixType<C>> Solver(HMatrix_, Eigen::ComputeEigenvectors);
        HMatrix_ = Solver.eigenvectors();
        Eigenvalues = Solver.eigenvalues(); // eigenvectors are ready
    }
}

template <bool C> MatrixType<C> const& HamiltonianPart::getMatrix() const {
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<MatrixType<C> const>(HMatrix);
}
template MatrixType<true> const& HamiltonianPart::getMatrix<true>() const;
template MatrixType<false> const& HamiltonianPart::getMatrix<false>() const;

template <bool C> MatrixType<C>& HamiltonianPart::getMatrix() {
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<MatrixType<C>>(HMatrix);
}
template MatrixType<true>& HamiltonianPart::getMatrix<true>();
template MatrixType<false>& HamiltonianPart::getMatrix<false>();

void HamiltonianPart::checkComputed() const {
    if(getStatus() < Computed)
        throw StatusMismatch("HamiltonianPart is not computed yet.");
}

RealType HamiltonianPart::getEigenValue(InnerQuantumState state) const {
    checkComputed();
    return Eigenvalues(static_cast<Eigen::Index>(state));
}

RealVectorType const& HamiltonianPart::getEigenValues() const {
    checkComputed();
    return Eigenvalues;
}

InnerQuantumState HamiltonianPart::getSize() const {
    return S.getBlockSize(Block);
}

template <bool C> VectorType<C> HamiltonianPart::getEigenState(InnerQuantumState state) const {
    checkComputed();
    return getMatrix<C>().col(state);
}
template VectorType<true> HamiltonianPart::getEigenState<true>(InnerQuantumState) const;
template VectorType<false> HamiltonianPart::getEigenState<false>(InnerQuantumState) const;

RealType HamiltonianPart::getMinimumEigenvalue() const {
    checkComputed();
    return Eigenvalues.minCoeff();
}

bool HamiltonianPart::reduce(RealType Cutoff) {
    checkComputed();

    Eigen::Index counter = 0;
    for(counter = 0; counter < Eigenvalues.size() && Eigenvalues[counter] <= Cutoff; ++counter)
        ;
    INFO("Left " << counter << " eigenvalues : ");

    if(counter) {
        INFO(Eigenvalues.head(counter) << "\n_________");
        Eigenvalues = Eigenvalues.head(counter);
        if(isComplex()) {
            auto& HMatrix_ = getMatrix<true>();
            HMatrix_ = HMatrix_.topLeftCorner(counter, counter);
        } else {
            auto& HMatrix_ = getMatrix<false>();
            HMatrix_ = HMatrix_.topLeftCorner(counter, counter);
        }
        return true;
    } else
        return false;
}

} // namespace Pomerol
