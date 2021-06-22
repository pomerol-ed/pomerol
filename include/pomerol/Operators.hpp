/** \file Operators.h
**  \brief Expressions with quantum-mechanical operators and functions to
*   construct them.
**
**  \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_OPERATORS_H
#define POMEROL_INCLUDE_OPERATORS_H

#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>
#include <libcommute/expression/hc.hpp>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Pomerol {
namespace Operators {

using libcommute::expression;
using libcommute::hc;

using libcommute::static_indices::c;
using libcommute::static_indices::c_dag;
using libcommute::static_indices::n;

namespace Detail {

//
// A C++11 compatible implementation of std::index_sequence
//
// Based on http://stackoverflow.com/a/17426611/410767 by Xeo
//

template <size_t... Ints> struct index_sequence {
    using type = index_sequence;
    using value_type = size_t;
    static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
};

// --------------------------------------------------------------
template <class Sequence1, class Sequence2> struct _merge_and_renumber;

template <size_t... I1, size_t... I2>
struct _merge_and_renumber<index_sequence<I1...>, index_sequence<I2...>>
    : index_sequence<I1..., (sizeof...(I1)+I2)...>
{};

// --------------------------------------------------------------

template <size_t N> struct make_index_sequence
    : _merge_and_renumber<typename make_index_sequence<N/2>::type,
                          typename make_index_sequence<N - N/2>::type>
{};

template<> struct make_index_sequence<0> : index_sequence<> {};
template<> struct make_index_sequence<1> : index_sequence<0> {};

//
// Apply a function to a tuple of arguments
//

template<typename T>
make_index_sequence<std::tuple_size<typename std::decay<T>::type>::value>
make_seq() { return {}; }

template<size_t N, typename T>
using element_t = typename std::tuple_element<N, typename std::decay<T>::type>::type;

template<typename F, typename ArgsT, size_t... Is>
auto apply_impl(F && f, ArgsT && args, index_sequence<Is...>) ->
    decltype(f(static_cast<element_t<Is, ArgsT>>(std::get<Is>(std::forward<ArgsT>(args)))...)) {
    return f(static_cast<element_t<Is, ArgsT>>(std::get<Is>(std::forward<ArgsT>(args)))...);
}

template<typename F, typename ArgsT>
auto apply(F && f, ArgsT && args) ->
    decltype(apply_impl(std::forward<F>(f), std::forward<ArgsT>(args), make_seq<ArgsT>())) {
    return apply_impl(std::forward<F>(f), std::forward<ArgsT>(args), make_seq<ArgsT>());
}

} // namespace Pomerol::Operators::detail

//
// Operator presets
//

template<typename... IndexTypes>
expression<double, IndexTypes...>
N(const std::vector<std::tuple<IndexTypes...>> &Indices) {
    expression<double, IndexTypes...> res;
    for(auto const& i : Indices)
        res += Detail::apply(n<double, IndexTypes...>, i);
    return res;
}

template<typename... IndexTypes>
expression<double, IndexTypes...>
Sz(const std::vector<std::tuple<IndexTypes...>> &SpinUpIndices,
   const std::vector<std::tuple<IndexTypes...>> &SpinDownIndices) {
    expression<double, IndexTypes...> res;
    for(auto const& i : SpinUpIndices)
        res += 0.5 * Detail::apply(n<double, IndexTypes...>, i);
    for(auto const& i : SpinDownIndices)
        res -= 0.5 * Detail::apply(n<double, IndexTypes...>, i);
    return res;
}

} // namespace Pomerol::Operators
} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_OPERATORS_H
