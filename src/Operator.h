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
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

namespace Pomerol{

/** A typedef for a generic single field operator. 
 * First bool = true, if it is a creation operator. Last argument is the index of the operator. */
typedef boost::tuple<bool, ParticleIndex> ElemOp;
/** Make an ElemOp and a vector<ElemOp> streamable. */
std::ostream& operator<< (std::ostream& output, const ElemOp& out);
std::ostream& operator<< (std::ostream& output, const std::vector<ElemOp>& out);

/** OpTerm is a sequence of elementary operators (std::vector) with a given matrix element. */
typedef boost::tuple<MelemType, std::vector<ElemOp> > OpTerm;
/** Make an OpTerm streamable. */
std::ostream& operator<< (std::ostream& output, const OpTerm& out);
/** Comparison operator. */
bool operator== (const OpTerm& lhs, const OpTerm& rhs);
/** A multiplication operator generates a term with a multiplication of matrix elements and 
 * a combined sequence of elementary operators, with first being the left hand side term. */
OpTerm operator*(const OpTerm& lhs, const OpTerm &rhs);
OpTerm operator*(const MelemType& lhs, const OpTerm &rhs);
OpTerm operator*(const OpTerm& lhs, const MelemType &rhs);

/** Operator represents a fermionic operator which is stored as a list of OpTerm's.
    This class is intended to store all operations, which are independent of the basis. */
class Operator
{
public:
    /** Checks the term for the consistency, e.g. Pauli principle and non-zero matrix element. */
    static bool checkTerm(const OpTerm& in);
protected:
    /** A set of Terms in the Operator. */
    std::list<OpTerm> Terms; // This will be inherited and used by classes

    /** Makes a swap of two adjacent operators in the term taking into account the anticommutation relation.
     * If operators anticommute, then a new term is generated and returned.
     * \param[in] position1 A position of the first operator to swap.
     * \param[in] position2 A position of the second operator to swap.
     * \param[in] force_ignore_commutation This forces to ignore all commutation relations and just to swap two operators and change the sign.
     * \param[out] A pair, which first argument is a result of the swap and the second is an operator, which
     * contains additional term, that may be produced while swapping.
     */
    static std::pair<OpTerm,Operator> elementary_swap_adjacent(const OpTerm &in, unsigned int position, bool force_ignore_commutation = false);
    /** A swap of two elements in the term. Has the same meaning as Operator::elementary_swap_adjacent, 
     * but can change any two operators in the term. 
     */
    static std::pair<OpTerm,Operator> elementary_swap(const OpTerm &in, unsigned int position1, unsigned position2, bool force_ignore_commutation = false);
    
    /** Returns a result of acting on a state to the right of an OpTerm
     * \param[in] ket A state to act on. 
     * \param[out] A pair of the resulting state and matrix element.
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
    Operator(const std::list<OpTerm> &Terms);
    /** Print all of the Terms. */
    void printAllTerms() const;
    /** Returns all Terms. */
    const std::list<OpTerm>& getTerms() const;

    /** Adds a term to the Operator. Done by checking (checkTerm) and push_back to the Terms list. 
      * \param[in] rhs A term to add.
      */
    Operator& operator+= (const OpTerm &rhs);
    /** Adds a term to the Operator with a matrix element multiplied by -1. */
    Operator& operator-= (const OpTerm &rhs);
    /** Adds all terms from the Operator rhs to the current Operator. */
    Operator& operator+= (const Operator &rhs);
    /** Adds all terms from the Operator rhs with the flipped sign of their matrix elements to the current Operator. */
    Operator& operator-= (const Operator &rhs);
    /** Returns a sum of current and rhs operator. */
    const Operator operator+(const Operator &rhs) const;
    /** Returns a sum of current and rhs operator with all matrix elements having a sign flipped. */
    const Operator operator-(const Operator &rhs) const;
    
    /** Returns an Operator, which has a list of Terms constructed by multiplying all of the terms of the current operator
     * to the ones in rhs. 
     */
    Operator operator*= (const Operator &rhs);
    /** Same as operator*=, but doesn't affect current Operator. */
    Operator operator* (const Operator &rhs) const;
    /** Multiplies all OpTerms of the operator by rhs. */
    Operator operator*= (const MelemType &rhs);
    /** Same as operator*=, but doesn't affect current Operator. */
    Operator operator* (const MelemType &rhs) const;
    /** Checks that two terms are equal. This is done, by rearranging all OpTerms in the Operators to the normal order,
      * sorting and comparing the resulting OpTerms. 
      */
    bool operator==(const Operator &rhs);
    /** Returns true if there are no OpTerms in this operator. */
    bool isEmpty() const;

    /** Returns an operator that contain the result of rearranging the OpTerms in this Operator according to a given rule.
     * \param[in] f A boost::function, that takes an OpTerm and returns the desired OpTerm.
     * \param[out] A resulting Operator. The original Operator is not changed. 
     */
    Operator rearrange(boost::function<std::vector<ElemOp> ( const std::vector<ElemOp> &in_f )> f ) const;
    /** Makes all Terms in the operator normal-ordered. */
    Operator getNormalOrdered() const;
    /** Reduces all terms with the same indices and order. */
    void reduce();
    /** Removes all terms from the list of terms with a given precision.
     *  \param[in] Terms A pointer to the list of pointers to Terms.
     */
    void prune(const RealType &Precision = std::numeric_limits<RealType>::epsilon());
    /** Sorts the list of the Term for comparison. */
    void sortTerms();
    
    /** Returns a matrix element of the operator.
     * \param[in] bra A state to the left of the operator.
     * \param[in] ket A state to the right of the operator.     
     * \param[out] Resulting matrix element.
     */
    virtual MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;

    /** Returns a result of acting of an operator on a state to the right of the operator.
     * \param[in] ket A state to act on.
     * \param[out] A map of states and corresponding matrix elements, which are the result of an action.
     */
    virtual std::map<FockState, MelemType> actRight(const FockState &ket) const;

    /** Returns an operator that is a commutator of the current operator and another one
     * \param[in] rhs An operator to calculate a commutator with.
     * \param[out] Resulting operator. 
     */
    Operator getCommutator(const Operator &rhs) const;
    
    /** Returns an operator that is an anticommutator of the current operator and another one
     * \param[in] rhs An operator to calculate an anticommutator with.
     * \param[out] Resulting operator. 
     */
    Operator getAntiCommutator(const Operator &rhs) const;

    /** Checks if current operator commutes with a given one. 
     * \param[in] rhs An operator to calculate a commutator with.
    */
    bool commutes(const Operator &rhs) const;
    
    /** Returns the total amount of Terms in the Operator. */
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
};

/** A small routine that returns bool, if the ElemOp contains a creation operator. Needed for STL algorithms. */
inline bool __isCdag(const ElemOp &in) { return in.get<0>() == 1; }
/** A comparison routine between two ElemOp's. */
inline bool __descendIndex(const ElemOp &in1, const ElemOp &in2) { return in1.get<1>()<in2.get<1>(); };

/** A routine that returns the normal ordered sequence of the ElemOp's in the current OpTerm. */
inline std::vector<ElemOp> NORMAL_ORDER ( const std::vector<ElemOp> &in_f)
{
    std::vector<ElemOp> out(in_f);
    std::vector<ElemOp>::iterator bound = std::partition(out.begin(), out.end(), __isCdag); 
    unsigned int index = std::distance(out.begin(), bound);
    std::sort(out.begin(), bound, __descendIndex); 
    std::sort(out.begin()+index, out.end(), __descendIndex); 
    return out;
}

/** A routine that returns the c^+c sequence of the ElemOp's in the current OpTerm. */
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

/** Generates a hash for an ElemOp. */
inline std::size_t ElemOphash_value(const ElemOp& in)
{
    std::size_t seed = 0;
    boost::hash_combine( seed, in.get<0>() );
    boost::hash_combine( seed, in.get<1>() );
    return seed;
}

#ifdef POMEROL_COMPLEX_MATRIX_ELEMENS
/** Comparison routines for two ElemOp. */
inline bool operator== (const ElemOp& lhs, const ElemOp& rhs){return (lhs.get<0>() == rhs.get<0>() && lhs.get<1>() == rhs.get<1>() );};
inline bool operator!= (const ElemOp& lhs, const ElemOp& rhs){return !(lhs==rhs);};
inline bool operator< (const ElemOp& lhs, const ElemOp& rhs)
{ 
    if (lhs.get<0>() == rhs.get<0>()) return (lhs.get<1>() < rhs.get<1>());
    return (lhs.get<0>() < rhs.get<0>());
};

inline bool operator> (const ElemOp& lhs, const ElemOp& rhs)
{
    return (lhs!=rhs && !(lhs<rhs));
}

inline bool operator< (const std::vector<ElemOp>& lhs, const std::vector<ElemOp>& rhs)
{
    if (lhs.size() < rhs.size()) return true;
    if (lhs.size() > rhs.size()) return false;
    for (int i=0; i<lhs.size(); ++i) if (lhs[i]<rhs[i]) return true; else if (lhs[i]>rhs[i]) return false;
    return false;
};

inline bool operator== (const std::vector<ElemOp>& lhs, const std::vector<ElemOp>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (int i=0; i<lhs.size(); ++i) if (lhs[i]!=rhs[i]) return false;
    return true;
};

/** A comparison routine for two OpTerms for complex valued matrix elements. */
inline
bool __compareOpTerms(const OpTerm& lhs, const OpTerm&rhs)
{
    const std::vector<ElemOp>& lv = lhs.get<1>();
    const std::vector<ElemOp>& rv = rhs.get<1>();
    if (lv == rv) return (std::abs(lhs.get<0>()) < std::abs(rhs.get<0>()));
    else return (lv < rv);
    //if (isEqual(lv,rv)) return (std::abs(lhs.get<0>()) < std::abs(rhs.get<0>()));
    //else return isLesser(lv, rv);
}
#endif


}; // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_OPERATOR_H
