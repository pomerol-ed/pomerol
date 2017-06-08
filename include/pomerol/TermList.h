//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2017 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2017 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.

/** \file include/TermList.h
** \brief List of terms in Lehmann representation
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef __INCLUDE_TERMLIST_H
#define __INCLUDE_TERMLIST_H

#include <set>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/serialization/set.hpp>

#include "Misc.h"

namespace Pomerol {

/** Container storing a list of terms (instances of TermType).
 *
 * The terms must allow sorting using a comparison function TermType::Compare.
 * Like terms (equivalent w.r.t. TermType::Compare) are automatically collected
 * and reduced to one term using operator+=().
 * A term T is considered negligible and is automatically removed from the
 * container if TermType::IsNegligible(T, current_number_of_terms + 1) evaluates to true.
 */
template<typename TermType> class TermList {

    /** Comparison function */
    typedef typename TermType::Compare Compare;
    /** Predicate for determining negligible terms */
    typedef typename TermType::IsNegligible IsNegligible;

    std::set<TermType, Compare> data;
    IsNegligible is_negligible;

public:

    /** Constructor.
     * \param[in] compare Compare predicate for the underlying std::set object
     * \param[in] is_negligible Predicate that determines whether a term can be neglected
     */
    TermList(Compare const& compare, IsNegligible const& is_negligible) :
        data(compare), is_negligible(is_negligible) {}

    /** Add a new term to the container
     * \param[in] term Term to be added
     */
    void add_term(TermType const& term) {
        typename std::set<TermType, Compare>::iterator it = data.find(term);
        if(it == data.end()) { // new term
            data.insert(term);
        } else {               // similar term
            TermType sum = *it;
            sum += term;
            data.erase(*it);
            if(!is_negligible(sum, data.size() + 1))
                data.insert(sum);
        }
    }

    /** Number of terms in the container */
    std::size_t size() const { return data.size(); }

    /** Remove all terms from the container */
    void clear() { data.clear(); }

    // Some pre-C++11 ugliness ...
#define MAKE_CALL_OPERATOR(N)                                               \
    template<BOOST_PP_ENUM_PARAMS(N, typename Arg)>                         \
    ComplexType operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, Arg, arg)) const {\
        ComplexType res = 0;                                                \
        for(typename std::set<TermType>::const_iterator it = data.begin();  \
            it != data.end(); ++it) {                                       \
            res += (*it)(BOOST_PP_ENUM_PARAMS(N, arg));                     \
        }                                                                   \
        return res;                                                         \
    }

    /** Pass arguments to operator() of each term in the container
     *  and return a sum of their return values.
     */
    MAKE_CALL_OPERATOR(1)
    MAKE_CALL_OPERATOR(2)
    MAKE_CALL_OPERATOR(3)
    MAKE_CALL_OPERATOR(4)
#undef MAKE_CALL_OPERATOR

    /** Boost.Serialization interface */
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        ar & data; ar & is_negligible;
    }

    /** Check that all terms in the container are properly ordered and are not negligible */
    bool check_terms() {
        if(size() == 0) return true;
        typename std::set<TermType>::const_iterator prev_it = data.begin();
        if(is_negligible(*prev_it, data.size() + 1)) return false;
        if(size() == 1) return true;

        Compare const& compare = data.key_comp();

        typename std::set<TermType>::const_iterator it = prev_it;
        for(++it; it != data.end(); ++it, ++prev_it) {
            if(is_negligible(*it, data.size() + 1) || !compare(*prev_it, *it))
                return false;
        }
        return true;
    }
};

}
#endif // endif :: #ifndef __INCLUDE_TERMLIST_H
