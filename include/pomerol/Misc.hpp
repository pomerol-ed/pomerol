//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Misc.hpp
/// \brief Declarations of the most basic types and macros.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifndef POMEROL_INCLUDE_POMEROL_MISC_HPP
#define POMEROL_INCLUDE_POMEROL_MISC_HPP

#include <pomerol/Version.hpp>

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/state_vector.hpp>

#ifndef DOXYGEN_SKIP
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif
#include <Eigen/Core>
#include <Eigen/Sparse>

#ifdef POMEROL_USE_OPENMP
#include <omp.h>
#endif

#if __cplusplus >= 201402L
#define POMEROL_DEPRECATED [[deprecated]]
#else
#define POMEROL_DEPRECATED
#endif

#include <array>
#include <complex>
#include <cstddef>
#include <iostream>
#include <type_traits>

/// The main namespace of the library.
namespace Pomerol {

/// \defgroup Basic Basic declarations
///@{

#ifndef DOXYGEN_SKIP
#define MSG_PREFIX __FILE__ << ":" << __LINE__ << ": "
#endif
#ifndef NDEBUG
/// Print a debugging message to the standard output with a source file name and line number annotation.
#define DEBUG(MSG) std::cout << MSG_PREFIX << MSG << '\n'
#else
#define DEBUG(MSG)
#endif
/// Print a message to the standard output.
#define INFO(MSG) std::cout << MSG << '\n'
/// Print a message without a trailing new line character to the standard output.
#define INFO_NONEWLINE(MSG) std::cout << MSG
/// Print a message to the standard error stream.
#define ERROR(MSG) std::cerr << MSG_PREFIX << MSG << '\n'

/// Real floating point type.
using RealType = double;
/// Complex floating point type.
using ComplexType = std::complex<double>;

/// Index of a single particle degree of freedom.
using ParticleIndex = unsigned int;

/// Index of a many-body state.
using QuantumState = libcommute::sv_index_type;

/// Dense complex matrix.
using ComplexMatrixType =
    Eigen::Matrix<ComplexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor>;

/// Element type of a real or complex matrix.
/// \tparam Complex Whether the matrix in question is complex.
template <bool Complex> using MelemType = typename std::conditional<Complex, ComplexType, RealType>::type;

/// Linear operator with a given type of coefficients.
/// \tparam ScalarType Coefficient  type.
template <typename ScalarType>
using LOperatorType = libcommute::loperator<ScalarType, libcommute::fermion, libcommute::boson>;

/// Linear operator with either real or complex coefficients.
/// \tparam Complex Whether the operator in question has complex coefficients.
template <bool Complex>
using LOperatorTypeRC = libcommute::loperator<MelemType<Complex>, libcommute::fermion, libcommute::boson>;

/// Dense real or complex matrix.
/// \tparam Complex Whether the matrix in question is complex.
template <bool Complex>
using MatrixType =
    Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor>;

/// Dense complex vector.
using ComplexVectorType = Eigen::Matrix<ComplexType, Eigen::Dynamic, 1, Eigen::AutoAlign>;
/// Dense real vector.
using RealVectorType = Eigen::Matrix<RealType, Eigen::Dynamic, 1, Eigen::AutoAlign>;

/// Dense real or complex vector.
/// \tparam Complex Whether the vector in question is complex.
template <bool Complex> using VectorType = Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic, 1, Eigen::AutoAlign>;

/// Sparse real or complex matrix with column-major storage.
/// \tparam Complex Whether the matrix in question is complex.
template <bool Complex> using ColMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::ColMajor>;
/// Sparse real or complex matrix with row-major storage.
/// \tparam Complex Whether the matrix in question is complex.
template <bool Complex> using RowMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::RowMajor>;

/// Imaginary unit \f$i\f$.
static ComplexType const I = ComplexType(0.0, 1.0); // 'static' to prevent linking problems

/// Permutation of 3 elements
struct Permutation3 {
    /// A permuted list of integers (0, 1, 2)
    std::array<std::size_t, 3> const perm;
    /// Signature of the permutation
    int const sign;
    bool operator==(Permutation3 const& rhs) const;
    bool operator!=(Permutation3 const& rhs) const;

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] p Permutation to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, Permutation3 const& p);
};
/// An array of all 3! = 6 permutations of 3 elements
extern std::array<Permutation3, 6> const permutations3;

/// Permutation of 4 elements
struct Permutation4 {
    /// A permuted list of integers (0, 1, 2, 3)
    std::array<std::size_t, 4> const perm;
    /// Signature of the permutation
    int const sign;
    bool operator==(Permutation4 const& rhs) const;
    bool operator!=(Permutation4 const& rhs) const;

    /// Output stream insertion operator.
    /// \param[out] os Output stream.
    /// \param[in] p Permutation to be inserted.
    /// \return Reference to the output stream.
    friend std::ostream& operator<<(std::ostream& os, Permutation4 const& p);
};
/// An array of all 4! = 24 permutations of 4 elements
extern std::array<Permutation4, 24> const permutations4;

/// Channel, in which a susceptibility function is defined.
enum Channel : int {
    PP, ///< Particle-particle channel.
    PH, ///< Particle-hole channel.
    xPH ///< Crossed particle-hole channel.
};
std::ostream& operator<<(std::ostream& os, Channel channel);

/// Hash function for real numbers that gives the same hash value for all
/// numbers falling into the same small interval (bin).
///
/// \param[in] x Value to be hashed.
/// \param[in] bin_size Size of the interval.
/// \return Hash value
std::size_t hash_binned_real(double x, double bin_size);

///@}

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_MISC_HPP
