#ifndef __INCLUDE_LIBCOMMUTEEIGEN_H
#define __INCLUDE_LIBCOMMUTEEIGEN_H

#include <libcommute/scalar_traits.hpp>
#include <libcommute/loperator/state_vector.hpp>

#include <Eigen/Core>

#include <type_traits>

namespace libcommute {

template<typename ScalarType,
         int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct element_type<Eigen::Matrix<ScalarType, Rows, Cols, Options, MaxRows, MaxCols>> {
    using type = ScalarType;
};

template<typename XprType, int BlockRows, int BlockCols, bool InnerPanel>
struct element_type<Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>> {
    using type = typename XprType::Scalar;
};

template<typename Derived>
inline auto get_element(Eigen::DenseBase<Derived> const& sv, sv_index_type n) ->
  typename Eigen::DenseBase<Derived>::Scalar const&
{
    return sv(n);
}

template<typename Derived, typename T>
inline void update_add_element(Eigen::DenseBase<Derived> & sv,
                               sv_index_type n,
                               T value) {
    sv(n) += value;
}

template<typename Derived>
inline void set_zeros(Eigen::DenseBase<Derived> & sv) {
    sv.setZero();
}

template<typename Derived, typename Functor>
inline void foreach(Eigen::DenseBase<Derived> const& sv, Functor&& f) {
    using ScalarType = typename Eigen::DenseBase<Derived>::Scalar;
    sv_index_type size = sv.size();
    for(sv_index_type n = 0; n < size; ++n) {
        auto const& a = sv(n);
        if(scalar_traits<ScalarType>::is_zero(a))
            continue;
        else
            f(n, a);
    }
}

} // end of namespace libcommute

#endif // endif :: #ifndef __INCLUDE_STATESCLASSIFICATION_H
