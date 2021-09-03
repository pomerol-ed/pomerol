//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file Misc.h
**    \brief Declares very common type names and macros.
**
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/
#ifndef POMEROL_INCLUDE_POMEROL_MISC_H
#define POMEROL_INCLUDE_POMEROL_MISC_H

#include <pomerol/Version.hpp>

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/state_vector.hpp>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Sparse>

#ifdef POMEROL_USE_OPENMP
#include <omp.h>
#endif

#include <array>
#include <complex>
#include <cstddef>
#include <iostream>
#include <type_traits>

namespace Pomerol {

#define MSG_PREFIX __FILE__ << ":" << __LINE__ << ": "
#ifndef NDEBUG
#define DEBUG(MSG) std::cout << MSG_PREFIX << MSG << std::endl
#else
#define DEBUG(MSG) NULL;
#endif
#define INFO(MSG) std::cout << MSG << std::endl
#define INFO_NONEWLINE(MSG) std::cout << MSG << std::flush
#define ERROR(MSG) std::cerr << MSG_PREFIX << MSG << std::endl

/** Real floating point type. */
using RealType = double;
/** Complex type. */
using ComplexType = std::complex<double>;

/** Index represents a combination of spin, orbital, and lattice indices **/
using ParticleIndex = unsigned int;

/** Each Quantum State in the finite system is associated with a number.
 * This works for any basis, including Fock and Hamiltonian eigenbasis.
 * The Fock States are converted naturally from bitsets to ints.
 **/
using QuantumState = libcommute::sv_index_type;

/** Index represents a combination of spin, orbital, and lattice indices **/
using ParticleIndex = unsigned int;

/** Dense complex matrix. */
using ComplexMatrixType =
    Eigen::Matrix<ComplexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor>;

template <bool Complex> using MelemType = typename std::conditional<Complex, ComplexType, RealType>::type;

template <typename ScalarType>
using LOperatorType = libcommute::loperator<ScalarType, libcommute::fermion, libcommute::boson>;

template <bool Complex>
using LOperatorTypeRC = libcommute::loperator<MelemType<Complex>, libcommute::fermion, libcommute::boson>;

template <bool Complex>
using MatrixType =
    Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor>;

/** Dense complex vector. */
using ComplexVectorType = Eigen::Matrix<ComplexType, Eigen::Dynamic, 1, Eigen::AutoAlign>;
/** Dense real vector. */
using RealVectorType = Eigen::Matrix<RealType, Eigen::Dynamic, 1, Eigen::AutoAlign>;

template <bool Complex> using VectorType = Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic, 1, Eigen::AutoAlign>;

/** Sparse complex matrix */
template <bool Complex> using ColMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::ColMajor>;
template <bool Complex> using RowMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::RowMajor>;

/** A short name for imaginary unit. */
static ComplexType const I = ComplexType(0.0, 1.0); // 'static' to prevent linking problems

/** Permutation of 3 elements */
struct Permutation3 {
    std::array<std::size_t, 3> const perm;
    int const sign;
    bool operator==(Permutation3 const& rhs) const;
    bool operator!=(Permutation3 const& rhs) const;
    friend std::ostream& operator<<(std::ostream& out, Permutation3 const& p);
};
extern std::array<Permutation3, 6> const permutations3;

/** Permutation of 4 elements */
struct Permutation4 {
    std::array<std::size_t, 4> const perm;
    int const sign;
    bool operator==(Permutation4 const& rhs) const;
    bool operator!=(Permutation4 const& rhs) const;
    friend std::ostream& operator<<(std::ostream& out, Permutation4 const& p);
};
extern std::array<Permutation4, 24> const permutations4;

} // namespace Pomerol

/**
 * \mainpage
 * The source code and fetch instructions are located at <a href="http://pomerol.googlecode.com">project's Google code page</a>.
 * \section   ref_API libpomerol API
 * The general sequence of a calculation is:
 * -    Define all indices of each mode of the system, i.e. site + spin indices ( by IndexClassification ).
 *      This is also a moment to discover symmetries against permutations of indices.
 * -    Define the Fock space - create QuantumStates and sort them into the blocks by their quantum numbers ( by StatesClassification ).
 * -    Enter blocks of the Hamiltonian ( by Hamiltonian and HamiltonianPart) and diagonalize them.
 * -    Find creation and annihilation operators in the eigenbasis of the Hamiltonian ( by FieldOperator and FieldOperatorPart ).
 * -    Calculate thermal quantities:
 *      -   The density matrix ( by DensityMatrix and DensityMatrixPart ).
 *      -   The Green's function ( by GreensFunction, GreensFunctionPart and GFContainer to store the values of GF for various index combinations ).
 *      -   The TwoParticle Greens Function ( by TwoParticleGF, TwoParticleGFPart and TwoParticleGFContainer ).
 * -    Calculcate the Vertex Function out of GF- and TwoParticleGF- containers - no work with Fock space is done ( by Vertex4 ).
 *
 * A hint: Refer to <a href="inherits.html">a Class Hierarchy</a> if provided
 *
 * \section ref_conventions Conventions
 *
 * \par Green's function
 * \f[
 *      G(\omega_n) = -\int_0^\beta \langle\mathbf{T}c_i(\tau)c^+_j(0)\rangle e^{i\omega_n\tau} d\tau
 * \f]
 *
 * \par Two-particle Green's function:
 * \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
 *   \int_0^\beta
 *     \langle\mathbf{T} c_i(\tau_1)c_j(\tau_2)c^+_k(\tau_3)c^+_l(0) \rangle
 *     \exp(i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3)
 *   d\tau_1 d\tau_2 d\tau_3
 * \f]
 *
 * \par The Wick part of a two-particle Green's function:
 * \f[
 * \chi^0_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) =
 * \beta\delta_{\omega_1\omega_4}\delta_{\omega_2\omega_3}G_{14}(\omega_1)G_{23}(\omega_2) -
 *  \beta\delta_{\omega_1\omega_3}\delta_{\omega_2\omega_4}G_{13}(\omega_1)G_{24}(\omega_2)
 * \f]
 * \par An irreducible vertex part:
 * \f[ \Gamma_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) \equiv
 *     \chi_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) -
 *     \chi^{0}_{1234}(\omega_1,\omega_2;\omega_3,\omega_4)
 * \f]
 *
 * \par An amputated irreducible vertex part:
 * \f[ \gamma_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) \equiv
 *     \sum_{1'2'3'4'}
 *     (G^{-1}(\omega_1))_{11'} (G^{-1}(\omega_2))_{22'}
 *     \Gamma_{1'2'3'4'}(\omega_1,\omega_2;\omega_3,\omega_4)
 *     (G^{-1}(\omega_3))_{3'3} (G^{-1}(\omega_4))_{4'4}
 * \f]
 *
 * \section Features
 * -    OpenMP support
 * \todo HDF5 storage ( read / write )
 * \section ref_authors Authors
 * -    Andrey Antipov <antipov[at]ct-qmc.org>
 * -    Igor Krivenko <igor[at]shg.ru>
 *
 * with a help from Mikhail Alejnikov, Alexey Rubtsov, Christoph Jung and Aljoscha Wilhelm
 *
 * We acknowledge <em>RRC Kurchatov Institute</em> for providing computing resources.
 *
 */

#endif // #ifndef POMEROL_INCLUDE_POMEROL_MISC_H
