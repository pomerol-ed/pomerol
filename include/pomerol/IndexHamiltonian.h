/** \file IndexHamiltonian.h
**  \brief Declaration of the IndexHamiltonian class - the Hamiltonian written in Index space.
**
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_INDEXHAMILTONIAN_H
#define __INCLUDE_INDEXHAMILTONIAN_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Lattice.h"
#include "Operator.h"

namespace Pomerol{

/* This class stores all matrix elements of a Hamiltonian in the index space. All terms have the ordering,
 * defined at the TERM_DEFAULT_SEQUENCE, which is by default taken as \f$ c^{\dagger} c c^{\dagger} c ... \f$. */
class IndexHamiltonian : public Operator
{
private:
    /** A pointer to the Lattice object. */
    const Lattice *L;
    /** A link to the IndexClassification object. */
    const IndexClassification &IndexInfo;
    /** A storage of Terms. Realized as a map of the order of the Term (number of operators) to the list of Terms. */
    //std::map <unsigned int, std::list<OpTerm*> > mapTerms;
public:
    /** Generates all Terms. */
    void prepare();
    /** Constructor. */
    IndexHamiltonian(const Lattice *L, const IndexClassification &Info);
    /** Gets all Terms of desired order. */
    //const std::list<OpTerm*> getTermsByOrder(unsigned int N) const;
    /** Print all Operator::Term s */
    //void printTerms(unsigned int order) const;
};

/** This function is used in the IndexHamiltonian prepare method.
 * This defines the default sequence of elements to be used in IndexHamiltonian OpTerm storage.. */
inline std::vector<bool>& TERM_DEFAULT_SEQUENCE(unsigned int N)
{
    static std::vector<bool> out;
    out.assign(N, false);
    for (unsigned int i=0; i<N; ++i) out[i]=(i+1)%2;
    return out;
}

}; // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_INDEXHAMILTONIAN_H
