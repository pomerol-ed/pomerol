//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file Operator.h
**  \brief Declarations of the Operator and Operator::Term classes.
** 
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_OPERATOR_H
#define __INCLUDE_OPERATOR_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Lattice.h"
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

namespace Pomerol{

/** This class represents an operator which is stored as a list of Terms */
class Operator
{
public:
    /** Declaration of a Term in the ParticleIndex space. */
    struct Term;
protected:
    /** A set of Terms in the Operator. */
    boost::shared_ptr<std::list<Operator::Term*> > Terms; // This will be inherited and used by classes
public:
    /** Empty constructor. */
    Operator();
    /** Constructor from the list of Terms. */
    Operator(boost::shared_ptr<std::list<Operator::Term*> > Terms);
    /** Print all of the Terms. */
    void printAllTerms() const;
    /** Returns all Terms. */
    boost::shared_ptr<std::list<Operator::Term*> > getTerms() const;

    /** Makes all Terms in the operator normal-ordered. */
    void makeNormalOrder();
    /** Reduces all terms with the same indices and order. */
    void reduce();
    /** Removes all terms from the list of terms with a given precision. */
    void prune(const RealType &Precision = std::numeric_limits<RealType>::epsilon());


    /** Returns a matrix element of the operator. */
    virtual MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;

    /** Returns a result of acting of an operator on a state
     * \param[in] ket A state to act on.
     * \param[out] A list of pairs of states and corresponding matrix elements, which are the result of an action.
     */
    virtual std::map<FockState, MelemType> actRight(const FockState &ket) const;

    /** Returns an operator that is a commutator of current operator and another one
     * \param[in] rhs An operator to calculate a commutator with.
     * \param[out] Resulting operator. */
    Operator getCommutator(const Operator &rhs) const;

    /** Checks if current operator commutes with a given one. 
     * \param[in] rhs An operator to calculate a commutator with.
    */
    bool commutes(const Operator &rhs) const;
    
    /** Destructor. */
    virtual ~Operator();
    friend std::ostream& operator<< (std::ostream& output, const Operator& out);
};

/** The Term in the ParticleIndex space is the same as the Lattice::Term, apart that it can be rearranged to the predefined sequence of operators
 * by using fermionic commutation relations. */
struct Operator::Term
{
friend class Operator;
protected:
    /** Number of operators in term. */
    const unsigned int N;
    /** Sequence of creation and annihilation operators. */
    std::vector<bool> OperatorSequence; 
    /** Array of ParticleIndices. */
    std::vector<ParticleIndex> Indices;
    /** Matrix element of Term. */
    MelemType Value;
private:

    /** Makes a swap of two adjacent operators in the term taking into account the anticommutation relation.
     * If operators anticommute, then a new term without these operators is returned.
     * \param[in] position A position of the first operator to swap
     * \param[in] force_ignore_commutation This forces to ignore all commutation relations and just to swap two operators and change the sign.
     * \param[out] Terms produced while swapping.
     */
    boost::shared_ptr<std::list<Operator::Term*> > elementary_swap(unsigned int position, bool force_ignore_commutation = false);

    /** Check if the Term is exactly the same as the other Term. */
    bool isExactlyEqual (const Term &rhs) const;

    /** Orders the indices in ascending/descending order. Works only for the term with the normal-ordered sequence of operators. */
    void reorder(bool ascend=true);
public:
    /** Reduces all terms with the same indices and order.
     *  \param[in] Terms A pointer to the list of pointers to Terms.
     */
    static void reduce(boost::shared_ptr<std::list<Operator::Term*> > Terms);
    /** Removes all terms from the list of terms with a given precision.
     *  \param[in] Terms A pointer to the list of pointers to Terms.
     */
    static void prune(boost::shared_ptr<std::list<Operator::Term*> > Terms, const RealType &Precision = std::numeric_limits<RealType>::epsilon());

    /** Rearranges operators in the term to a desired sequence. Warning! This operation is not unique
     * \param[in] DesiredSequence A sequence of operators ( represented as a vector of bool ) to rearrange the term.
     */
    boost::shared_ptr<std::list<Operator::Term*> > rearrange(const std::vector<bool> & DesiredSequence); 
    /** Rearranges a term to the normal order (с^+ to the left, c to the right). Makes ascending order of the indices. */
    boost::shared_ptr<std::list<Operator::Term*> > makeNormalOrder();
    /** Return amount of operators in Term. */
    unsigned int getN();
    /** Returns a matrix element of the Operator::Term. */
    MelemType getMatrixElement(const FockState &bra, const FockState &ket);
    /** Check if the Term is the same as the other Term. */
    bool operator==(const Term &rhs) const;
    /** Check if the Term commutes with the other Term. */
    bool commutes(const Term &rhs) const;

    /** Calculate a commutator with another Operator::Term.
     * \param[in] rhs An operator to calculate the commutator.
     * \param[out] A resulting list of two terms. 
     */
    boost::shared_ptr<std::list<Operator::Term*> > getCommutator(const Operator::Term &rhs) const;
    /** Returns a result of acting on a state by a Term
     * \param[in] ket A state to act on. 
     * \param[out] A pair of Resulting state and matrix element.
     */
    boost::tuple<FockState, MelemType> actRight(const FockState &ket);

    /** Constructor
     * \param[in] N Total amount of operators in the term.
     * \param[in] Sequence Sequence of creation/annihilation operators in the term. True goes for creation, false - for annihilation.
     * \param[in] Indices Corresponding indices of the creation/annihilation operators.
     */
    Term (const unsigned int N, const std::vector<bool>&  Sequence, const std::vector<ParticleIndex>& Indices, MelemType Value);
    /** Exception - wrong operation with labels. */
    class exWrongLabel : public std::exception { virtual const char* what() const throw(); };
    /** Exception - wrong operation with bool sequence. */
    class exWrongOpSequence : public std::exception { virtual const char* what() const throw(); };
    /** Exception - Matrix element of term vanishes. */
    class exMelemVanishes : public std::exception { virtual const char* what() const throw(); };

/** Make the Term printable */
friend std::ostream& operator<< (std::ostream& output, const Term& out);
};

}; // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_OPERATOR_H
