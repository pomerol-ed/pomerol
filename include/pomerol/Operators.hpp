//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/pomerol/Operators.hpp
/// \brief Expressions with quantum-mechanical operators and functions to construct them.
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_OPERATORS_HPP
#define POMEROL_INCLUDE_OPERATORS_HPP

#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>
#include <libcommute/expression/hc.hpp>

#include <cstdlib>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Pomerol {
namespace Operators {

/// \defgroup Operators Operator expressions
///@{

/// \namespace Pomerol::Operators
/// \brief Expressions of quantum-mechanical operators.
///
/// The \ref Pomerol::Operators namespace imports a few types and functions from
/// <a href="https://krivenko.github.io/libcommute/">libcommute</a>. These include
/// \li <a href="https://krivenko.github.io/libcommute/expression/expression.html">
///     The polynomial expression object, libcommute::expression</a>.
/// \li <a href="https://krivenko.github.io/libcommute/expression/expression.html#pm-h-c-notation">
///     The plus/minus Hermitian conjugate placeholder, libcommute::hc</a>.
/// \li <a href="https://krivenko.github.io/libcommute/expression/factories.html#statically-typed-indices">
///     Factory functions for fermionic creation, annihilation and occupation operators
///     (\f$ \hat c^\dagger, \hat c, \hat n = \hat c^\dagger \hat c\f$) with statically typed indices</a>.
/// \li <a href="https://krivenko.github.io/libcommute/expression/factories.html#statically-typed-indices">
///     Factory functions for bosonic creation and annihilation operators
///     (\f$ \hat a^\dagger, \hat a\f$) with statically typed indices</a>.
///
/// There are also a few additional factory functions defined in this namespace.

using libcommute::expression;
using libcommute::hc;

using libcommute::static_indices::c;
using libcommute::static_indices::c_dag;
using libcommute::static_indices::n;
using libcommute::static_indices::a;
using libcommute::static_indices::a_dag;

namespace Detail {

//
// A C++11 compatible implementation of std::index_sequence
//
// Based on http://stackoverflow.com/a/17426611/410767 by Xeo
//

template <std::size_t... Ints> struct index_sequence {
    using type = index_sequence;
    using value_type = std::size_t;
    static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
};

// --------------------------------------------------------------
template <class Sequence1, class Sequence2> struct _merge_and_renumber;

template <std::size_t... I1, std::size_t... I2>
struct _merge_and_renumber<index_sequence<I1...>, index_sequence<I2...>>
    : index_sequence<I1..., (sizeof...(I1) + I2)...> {};

// --------------------------------------------------------------

template <std::size_t N>
struct make_index_sequence
    : _merge_and_renumber<typename make_index_sequence<N / 2>::type, typename make_index_sequence<N - N / 2>::type> {};

template <> struct make_index_sequence<0> : index_sequence<> {};
template <> struct make_index_sequence<1> : index_sequence<0> {};

//
// Apply a function to a tuple of arguments
//

template <typename T> make_index_sequence<std::tuple_size<typename std::decay<T>::type>::value> make_seq() {
    return {};
}

template <std::size_t N, typename T>
using element_t = typename std::tuple_element<N, typename std::decay<T>::type>::type;

template <typename F, typename ArgsT, std::size_t... Is>
auto apply_impl(F&& f, ArgsT&& args, index_sequence<Is...>)
    -> decltype(f(static_cast<element_t<Is, ArgsT>>(std::get<Is>(std::forward<ArgsT>(args)))...)) {
    return f(static_cast<element_t<Is, ArgsT>>(std::get<Is>(std::forward<ArgsT>(args)))...);
}

template <typename F, typename ArgsT>
auto apply(F&& f, ArgsT&& args)
    -> decltype(apply_impl(std::forward<F>(f), std::forward<ArgsT>(args), make_seq<ArgsT>())) {
    return apply_impl(std::forward<F>(f), std::forward<ArgsT>(args), make_seq<ArgsT>());
}

} // namespace Detail

//
// Operator presets
//

/// Construct a real-valued expression for the full occupation number operator
/// \f[
///   \hat N = \sum_{i\in\{I\}} \hat n_i.
/// \f]
/// \tparam IndexTypes Types of indices of creation/annihilation operators in the resulting expression.
/// \param[in] Indices List of index tuples corresponding to the selected degrees of freedom \f$\{I\}\f$.
/// \return Constructed expression.
template <typename... IndexTypes>
expression<double, IndexTypes...> N(std::vector<std::tuple<IndexTypes...>> const& Indices) {
    expression<double, IndexTypes...> res;
    for(auto const& i : Indices)
        res += Detail::apply(n<double, IndexTypes...>, i); // cppcheck-suppress useStlAlgorithm
    return res;
}

/// Construct a real-valued expression for the full spin z-projection operator
/// \f[
///   \hat S_z = \frac{1}{2}\sum_{i\in\{I_\uparrow\}}\hat n_i - \frac{1}{2}\sum_{i\in\{I_\downarrow\}}\hat n_i.
/// \f]
/// \tparam IndexTypes Types of indices of creation/annihilation operators in the resulting expression.
/// \param[in] SpinUpIndices List of index tuples corresponding to the spin-up degrees of freedom
/// \f$\{I_\uparrow\}\f$.
/// \param[in] SpinDownIndices List of index tuples corresponding to the spin-down degrees of freedom
/// \f$\{I_\downarrow\}\f$.
/// \return Constructed expression.
template <typename... IndexTypes>
expression<double, IndexTypes...> Sz(std::vector<std::tuple<IndexTypes...>> const& SpinUpIndices,
                                     std::vector<std::tuple<IndexTypes...>> const& SpinDownIndices) {
    expression<double, IndexTypes...> res;
    for(auto const& i : SpinUpIndices)
        res += 0.5 * Detail::apply(n<double, IndexTypes...>, i); // cppcheck-suppress useStlAlgorithm
    for(auto const& i : SpinDownIndices)
        res -= 0.5 * Detail::apply(n<double, IndexTypes...>, i); // cppcheck-suppress useStlAlgorithm
    return res;
}

///@}

} // namespace Operators
} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_OPERATORS_HPP
