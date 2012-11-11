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
**  \brief Declarations of the Operator and OpTerm classes.
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
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

namespace Pomerol{

/** A typedef for a generic single field operator. 
 * First bool = true, if it is a creation operator. Last argument is the index of the operator. */
typedef boost::tuple<bool, ParticleIndex> ElemOp;
std::ostream& operator<< (std::ostream& output, const ElemOp& out);
std::ostream& operator<< (std::ostream& output, const std::vector<ElemOp>& out);
/** A vector of elementary operators generate a term with a given matrix element. */
typedef boost::tuple<MelemType, std::vector<ElemOp> > OpTerm;
std::ostream& operator<< (std::ostream& output, const OpTerm& out);
bool operator== (const OpTerm& lhs, const OpTerm& rhs);
OpTerm operator*(const OpTerm& lhs, const OpTerm &rhs);

/** This class represents an operator which is stored as a list of Terms */
class Operator
{
public:
    /** Checks the term for the consistency. */
    static bool checkTerm(const OpTerm& in);
protected:
    /** A set of Terms in the Operator. */
    boost::shared_ptr<std::list<OpTerm> > Terms; // This will be inherited and used by classes

    /** Makes a swap of two adjacent operators in the term taking into account the anticommutation relation.
     * If operators anticommute, then a new term without these operators is returned.
     * \param[in] position1 A position of the first operator to swap.
     * \param[in] position2 A position of the second operator to swap.
     * \param[in] force_ignore_commutation This forces to ignore all commutation relations and just to swap two operators and change the sign.
     * \param[out] Terms produced while swapping.
     */
    static std::pair<OpTerm,Operator> elementary_swap(const OpTerm &in, unsigned int position1, unsigned position2, bool force_ignore_commutation = false);
    static std::pair<OpTerm,Operator> elementary_swap_adjacent(const OpTerm &in, unsigned int position, bool force_ignore_commutation = false);

    /** Returns a result of acting on a state by an OpTerm
     * \param[in] ket A state to act on. 
     * \param[out] A pair of Resulting state and matrix element.
     */
    static boost::tuple<FockState,MelemType> actRight(const OpTerm &in, const FockState &ket);
    //static boost::tuple<FockState,MelemType> actLeft(const OpTerm &in, const FockState &bra);


    /** Rearranges the term according to a function f, which takes the vector of operators and returns a desired vector of operators.
     *  Warning! This operation might not be unique.
     * \param[in] f A Boost::function that rearranges the operators in the term.
     * \param[in] in An OpTerm to rearrange.
     */


public:
    /** Empty constructor. */
    Operator();
    /** Constructor from the single term. */
    Operator(const OpTerm& term);
    /** Constructor from the list of Terms. */
    Operator(boost::shared_ptr<std::list<OpTerm> > Terms);
    /** Print all of the Terms. */
    void printAllTerms() const;
    /** Returns all Terms. */
    boost::shared_ptr<std::list<OpTerm> > getTerms() const;

    void add (const OpTerm &rhs);
    Operator& operator+= (const OpTerm &rhs);
    Operator& operator+= (const Operator &rhs);
    const Operator operator+(const Operator &rhs) const;
    Operator operator*= (const Operator &rhs);
    Operator operator* (const Operator &rhs) const;
    bool operator==(const Operator &rhs);
    bool isEmpty() const;

    Operator rearrange(boost::function<std::vector<ElemOp> ( const std::vector<ElemOp> &in_f )> f ) const;
    /** Makes all Terms in the operator normal-ordered. */
    Operator getNormalOrdered() const;
    /** Reduces all terms with the same indices and order. */
    void reduce();
    /** Removes all terms from the list of terms with a given precision.
     *  \param[in] Terms A pointer to the list of pointers to Terms.
     */
    void prune(const RealType &Precision = std::numeric_limits<RealType>::epsilon());
    /** Reduces all terms with the same indices and order.
     *  \param[in] Terms A pointer to the list of pointers to Terms.
     */

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
    
    const unsigned int getNTerms() const;
    /** Make the Operator printable */
    friend std::ostream& operator<< (std::ostream& output, const Operator& out);

    /** Exception - wrong operation with labels. */
    class exWrongLabel : public std::exception { virtual const char* what() const throw(); };
    /** Exception - wrong operation with bool sequence. */
    class exWrongOpSequence : public std::exception { virtual const char* what() const throw(); };
    /** Exception - Matrix element of term vanishes. */
    class exMelemVanishes : public std::exception { virtual const char* what() const throw(); };

    /** Destructor. */
    virtual ~Operator();
    friend std::ostream& operator<< (std::ostream& output, const Operator& out);
};


inline bool __isCdag(const ElemOp &in) { return in.get<0>() == 1; }
inline bool __descendIndex(const ElemOp &in1, const ElemOp &in2) { return in1.get<1>()<in2.get<1>(); };

inline std::vector<ElemOp> NORMAL_ORDER ( const std::vector<ElemOp> &in_f)
{
    std::vector<ElemOp> out(in_f);
    std::vector<ElemOp>::iterator bound = std::partition(out.begin(), out.end(), __isCdag); 
    unsigned int index = std::distance(out.begin(), bound);
    std::sort(out.begin(), bound, __descendIndex); 
    std::sort(out.begin()+index, out.end(), __descendIndex); 
    return out;
}

inline std::vector<ElemOp> CDAG_C ( const std::vector<ElemOp> &in_f)
{
    std::vector<ElemOp> out(in_f);
    unsigned int N = in_f.size();
    std::vector<ElemOp>::iterator bound = std::partition(out.begin(), out.end(), __isCdag); 
    unsigned int index = std::distance(out.begin(), bound);
    std::sort(out.begin(), bound, __descendIndex); 
    std::sort(out.begin()+index, out.end(), __descendIndex); 
    std::vector<ElemOp> out2(out);
    unsigned int j1=0, j2=index, out2_index=0;
    for (out2_index=0; out2_index<N; out2_index++) {
        if (j1<index && j2 < N) 
            if (out2_index%2==0) { out2[out2_index]=out[j1]; j1++; }
            else { out2[out2_index]=out[j2]; j2++; }
        else if (j1 < index) { out2[out2_index]=out[j1]; j1++; }
             else { out2[out2_index]=out[j2]; j2++; };
    }
    return out;
}


inline std::size_t ElemOphash_value(const ElemOp& in)
{
    std::size_t seed = 0;
    boost::hash_combine( seed, in.get<0>() );
    boost::hash_combine( seed, in.get<1>() );
    return seed;
}

}; // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_OPERATOR_H
