//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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


/** \file IndexHamiltonian.h
**  \brief Declaration of the IndexHamiltonian class - the Hamiltonian written in Index space.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_INDEXHAMILTONIAN_H
#define __INCLUDE_INDEXHAMILTONIAN_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Lattice.h"
#include "Operator.h"
#include <boost/shared_ptr.hpp>

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
    std::map <unsigned int, std::list<Operator::Term*> > mapTerms;
public:
    /** Generates all Terms. */
    void prepare();
    /** Constructor. */
    IndexHamiltonian(const Lattice *L, const IndexClassification &Info);
    /** Gets all Terms of desired order. */
    const std::list<Operator::Term*> getTermsByOrder(unsigned int N) const;
    /** Print all Operator::Term s */
    void printTerms(unsigned int order) const;
};

/** This function is used in the IndexHamiltonian prepare method. 
 * This defines the default sequence of elements to be used in IndexHamiltonian Operator::Term storage.. */
inline std::vector<bool>& TERM_DEFAULT_SEQUENCE(unsigned int N)
{
    static std::vector<bool> out;
    out.assign(N, false);
    for (unsigned int i=0; i<N; ++i) out[i]=(i+1)%2;
    return out;
}

}; // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_INDEXHAMILTONIAN_H
