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

namespace Pomerol{

/* This class stores all matrix elements of a Hamiltonian in the index space. All terms have the ordering, 
 * defined at the TERM_DEFAULT_SEQUENCE, which is by default taken as \f$ c^{\dagger} c c^{\dagger} c ... \f$. */
class IndexHamiltonian
{
public:
    /** Declaration of a Term in the ParticleIndex space. */
    struct Term;
private:
    /** A pointer to the Lattice object. */
    const Lattice *L;
    /** A link to the IndexClassification object. */
    const IndexClassification &IndexInfo;
    /** A storage of Terms. Realized as a map of the order of the Term (number of operators) to the list of Terms. */
    std::map <unsigned int, std::list<Term*> > Terms;
public:
    /** Generates all Terms. */
    void prepare();
    /** Constructor. */
    IndexHamiltonian(const Lattice *L, const IndexClassification &Info);
    /** Gets all Terms of desired order. */
    const std::list<Term*> getTerms(unsigned int N) const;
    /** Print all IndexHamiltonian::Term s */
    void printTerms(unsigned int order) const;
};

/** The Term in the ParticleIndex space is the same as the Lattice::Term, apart that it can be rearranged to the predefined sequence of operators
 * by using fermionic commutation relations. */
struct IndexHamiltonian::Term
{
friend class IndexHamiltonian;
protected:
    /** Number of operators in term. */
    const unsigned int Order;
    /** Sequence of creation and annihilation operators. */
    std::vector<bool> OperatorSequence; 
    /** Array of ParticleIndices. */
    std::vector<ParticleIndex> Indices;
    /** Matrix element of Term. */
    RealType Value;
private:

    /** Makes a swap of two adjacent operators in the term taking into account the anticommutation relation.
     * If operators anticommute, then a new term without these operators is returned.
     * \param[in] position A position of the first operator to swap
     * \param[in] force_ignore_commutation This forces to ignore all commutation relations and just to swap two operators and change the sign.
     * \param[out] Terms produced while swapping.
     */
    std::list<IndexHamiltonian::Term*> elementary_swap(unsigned int position, bool force_ignore_commutation = false);
public:
    /** Rearranges operators in the term to a desired sequence. 
     * \param[in] DesiredSequence A sequence of operators ( represented as a vector of bool ) to rearrange the term.
     */
    std::list<IndexHamiltonian::Term*> rearrange(const std::vector<bool> & DesiredSequence); 

    /** Constructor
     * \param[in] Order Total amount of operators in the term.
     * \param[in] Sequence Sequence of creation/annihilation operators in the term. True goes for creation, false - for annihilation.
     * \param[in] Indices Corresponding indices of the creation/annihilation operators.
     */
    Term (const unsigned int Order, const std::vector<bool>&  Sequence, const std::vector<ParticleIndex>& Indices, RealType Value);
    /** Exception - wrong operation with labels. */
    class exWrongLabel : public std::exception { virtual const char* what() const throw(); };
    /** Exception - wrong operation with bool sequence. */
    class exWrongOpSequence : public std::exception { virtual const char* what() const throw(); };

/** Make the Term printable */
friend std::ostream& operator<< (std::ostream& output, const Term& out);
};

/** This function is used in the IndexHamiltonian prepare method. 
 * This defines the default sequence of elements to be used in IndexHamiltonian IndexHamiltonian::Term storage.. */
inline std::vector<bool>& TERM_DEFAULT_SEQUENCE(unsigned int N)
{
    static std::vector<bool> out;
    out.assign(N, false);
    for (unsigned int i=0; i<N; ++i) out[i]=(i+1)%2;
    return out;
}

}; // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_INDEXHAMILTONIAN_H
