/** \file include/TermList.h
** \brief List of terms in Lehmann representation
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_TERMLIST_H
#define POMEROL_INCLUDE_TERMLIST_H

#include "mpi_dispatcher/misc.hpp"

#include "Misc.h"

#include <set>
#include <utility>
#include <vector>

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
        auto it = data.find(term);
        if(it == data.end()) { // new term
            data.emplace(term);
        } else {               // similar term
            TermType sum = *it;
            sum += term;
            data.erase(*it);
            if(!is_negligible(sum, data.size() + 1))
                data.emplace(sum);
        }
    }

    /** Number of terms in the container */
    std::size_t size() const { return data.size(); }

    /** Remove all terms from the container */
    void clear() { data.clear(); }

    /** Access set of terms */
    const std::set<TermType, Compare>& as_set() const { return data; }

    /** Access the is_negligible object */
    const IsNegligible& get_is_negligible() const { return is_negligible; }

    /** Pass arguments to operator() of each term in the container
     *  and return a sum of their return values.
     */
    template<typename... Args>
    ComplexType operator()(Args&&... args) const {
        ComplexType res = 0;
        for(auto const& t : data)
            res += t(std::forward<Args>(args)...);
        return res;
    }

    void broadcast(const MPI_Comm &comm, int root) {
        auto comp = data.key_comp();
        comp.broadcast(comm, root);

        long n_terms;
        if(pMPI::rank(comm) == root) { // Broadcast the terms from this process
            n_terms = data.size();
            std::vector<TermType> v(data.begin(), data.end());
            MPI_Bcast(&n_terms, 1, MPI_LONG, root, comm);
            MPI_Bcast(v.data(), v.size(), TermType::mpi_datatype(), root, comm);
        } else { // Receive terms
            MPI_Bcast(&n_terms, 1, MPI_LONG, root, comm);
            std::vector<TermType> v(n_terms);
            MPI_Bcast(v.data(), v.size(), TermType::mpi_datatype(), root, comm);
            data = std::set<TermType, Compare>(v.begin(), v.end(), comp);
        }

        is_negligible.broadcast(comm, root);
    }

    /** Check that all terms in the container are properly ordered and are not negligible */
    bool check_terms() {
        if(size() == 0) return true;
        auto prev_it = data.begin();
        if(is_negligible(*prev_it, data.size() + 1)) return false;
        if(size() == 1) return true;

        Compare const& compare = data.key_comp();

        auto it = prev_it;
        for(++it; it != data.end(); ++it, ++prev_it) {
            if(is_negligible(*it, data.size() + 1) || !compare(*prev_it, *it))
                return false;
        }
        return true;
    }
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TERMLIST_H
