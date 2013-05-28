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


/** \file Operator.cpp
**  \brief Implementation of the Operator, Operator::Term classes
** 
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "Operator.h"
#include <algorithm>
#include <iterator>
#include <boost/tuple/tuple.hpp>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/utility/swap.hpp"
#include <boost/functional/hash.hpp>

namespace Pomerol{

std::ostream& operator<< (std::ostream& output, const ElemOp& out)
{
    output << (out.get<0>()?"c^{+}":"c") << "_" << out.get<1>();
    return output;
}

std::ostream& operator<< (std::ostream& output, const std::vector<ElemOp>& out)
{
    for (unsigned int i=0; i<out.size(); ++i)  output << out[i]; 
    return output;
}

std::ostream& operator<< (std::ostream& output, const OpTerm& out)
{
    output << out.get<0>();
    const std::vector<ElemOp> &elems = out.get<1>();
    for (unsigned int i=0; i<elems.size(); ++i)  output << elems[i]; 
    return output;
}

bool operator== (const OpTerm& lhs, const OpTerm& rhs)
{
    return (std::abs(lhs.get<0>()-rhs.get<0>())<std::numeric_limits<RealType>::epsilon() && lhs.get<1>() == rhs.get<1>());
}

OpTerm operator*(const OpTerm& lhs, const OpTerm &rhs)
{
    OpTerm out(lhs);
    out.get<0>() = lhs.get<0>()*rhs.get<0>();
    out.get<1>().insert(out.get<1>().end(),rhs.get<1>().begin(),rhs.get<1>().end());
    return out;
}

OpTerm operator*(const OpTerm& lhs, const MelemType &rhs)
{
    OpTerm out(lhs);
    out.get<0>() = lhs.get<0>()*rhs;
    return out;
}

OpTerm operator*(const MelemType& lhs, const OpTerm &rhs)
{
    return rhs*lhs;
}

bool Operator::checkTerm(const OpTerm &in)
{
    if (std::abs(in.get<0>()) < std::numeric_limits<RealType>::epsilon()) return false;
    const std::vector<ElemOp> &ops = in.get<1>(); 
    unsigned int N=ops.size();
    for (unsigned int i=0; i<N; ++i) {  
        bool c; ParticleIndex ind;
        boost::tie(c,ind) = ops[i]; 
        int count_index=2*c-1; // This determines how many times current index Indices[i] is found. c^+ gives +1, c gives -1.
        for (unsigned int j=i+1; j<N; ++j) {
            bool c2; ParticleIndex ind2;
            boost::tie(c2,ind2)=ops[j];
            if (ind2 == ind ) {  //same indices
                count_index+=2*c2-1; 
                if ( count_index > 1 || count_index < -1 ) { 
                    // ERROR("This term " << in << " vanishes. "); 
                    return false; 
                    }; 
                } 
            }
        };    
    return true;
}

//
// Operator
//

Operator::Operator()
{
    //Terms.reset( new std::list<OpTerm> );
}

Operator::Operator(const std::list<OpTerm>& Terms2) : Terms(Terms2)
{
    Terms.remove_if(this->checkTerm);
}

Operator::Operator(const OpTerm& term)
{
    if (checkTerm(term)) Terms.push_back(term);
}

const std::list<OpTerm>& Operator::getTerms() const
{
    return Terms;
}

const unsigned int Operator::getNTerms() const
{
    return Terms.size();
}

Operator::~Operator()
{
    //~Terms();
}

bool Operator::isEmpty() const
{
    return (Terms.size()==0);
}

std::ostream& operator<< (std::ostream& output, const Operator& out)
{
    for (std::list<OpTerm>::const_iterator it = out.Terms.begin(); it!=out.Terms.end(); it++) {
        if (it!=out.Terms.begin())  output << " + ";
        output << *it;
        };
    return output;
}

void Operator::printAllTerms() const
{
    INFO(*this);
}

Operator& Operator::operator+=(const Operator &rhs)
{
    if ( rhs.Terms.size() ) Terms.insert(Terms.end(), rhs.Terms.begin(), rhs.Terms.end()); 
    return *this;
}

Operator& Operator::operator+=(const OpTerm &rhs)
{
    if (checkTerm(rhs)) Terms.push_back(rhs);
    return *this;
}

Operator& Operator::operator-=(const OpTerm &rhs)
{
    if (checkTerm(rhs)) Terms.push_back((-1.0)*rhs);
    return *this;
}

Operator& Operator::operator-=(const Operator &rhs)
{
    (*this)+=rhs*(-1.);
    return *this;
}

const Operator Operator::operator+(const Operator &rhs) const
{
    Operator out(*this);
    out+=rhs;
    return out;
}

const Operator Operator::operator-(const Operator &rhs) const
{
    Operator out(*this);
    out-=rhs;
    return out;
}

Operator Operator::operator*=(const Operator &rhs)
{
    Operator out;
    for ( std::list<OpTerm>::iterator term_lhs_it = Terms.begin(); term_lhs_it != Terms.end(); term_lhs_it++) {
        for ( std::list<OpTerm>::const_iterator term_rhs_it = rhs.Terms.begin(); term_rhs_it != rhs.Terms.end(); term_rhs_it++) {
            //DEBUG(*term_lhs_it);
            //DEBUG(*term_rhs_it);
            out+=((*term_lhs_it)*(*term_rhs_it));
        }
    }
    *this = out;
    return *this;
}

Operator Operator::operator*(const Operator &rhs) const
{
    Operator out=(*this);
    out*=rhs;
    return out;
}

Operator Operator::operator*=(const MelemType &rhs)
{
    for ( std::list<OpTerm>::iterator term_lhs_it = Terms.begin(); term_lhs_it != Terms.end(); term_lhs_it++) {
            *term_lhs_it=(*term_lhs_it)*rhs;
        }
    return *this;
}

Operator Operator::operator*(const MelemType &rhs) const
{
    Operator out=(*this);
    out*=rhs;
    return out;
}

std::pair<OpTerm, Operator> Operator::elementary_swap_adjacent(const OpTerm &in, unsigned int position, bool force_ignore_commutation)
{
    assert (position < in.get<1>().size()-1);
    Operator out;
    OpTerm out_term(in);
    std::vector<ElemOp>& ops = out_term.get<1>();  MelemType Value=out_term.get<0>();
    if ( ops[position].get<1>() != ops[position+1].get<1>() || force_ignore_commutation ) {
        out_term.get<0>()*=(-1.); 
        boost::swap(ops[position], ops[position+1]);
        }
    else {
        if (ops[position].get<0>() == ops[position+1].get<0>()) throw (exWrongOpSequence()); 
        OpTerm term;
        //DEBUG("Here!!!");
        term.get<1>().insert(term.get<1>().end(), ops.begin(), ops.begin()+position);
        term.get<1>().insert(term.get<1>().end(), ops.begin()+position+2, ops.end());
        term.get<0>() = Value;
        //DEBUG("Now making swap");
        out_term = elementary_swap_adjacent(out_term, position, true).first;
        //DEBUG("received " << out_term << " out of a swap at pos " << position << " of " << in);
        out+=term;
         };
    //DEBUG(in << "--->" << out_term);
    return std::make_pair(out_term,out);
}

std::pair<OpTerm,Operator> Operator::elementary_swap(const OpTerm &in, unsigned int position1, unsigned int position2, bool force_ignore_commutation)
{
    Operator out_op;
    OpTerm out_term(in);
    std::vector<ElemOp>& ops = out_term.get<1>();  MelemType Value=out_term.get<0>();
    if (position2 == position1) return std::make_pair(out_term,out_op);
    if (position2 < position1) std::swap(position2, position1);
    for (unsigned int i=position1; i<position2; i++) { 
        std::pair<OpTerm,Operator> out_pair = elementary_swap_adjacent(out_term, i);
        out_term = out_pair.first;
        out_op+=out_pair.second;
        }
    for (unsigned int i=position2-2; i>=position1; i++) {
        std::pair<OpTerm,Operator> out_pair = elementary_swap_adjacent(out_term, i);
        out_term = out_pair.first;
        out_op+=out_pair.second;
        }
    return std::make_pair(out_term,out_op);
}

Operator Operator::rearrange(boost::function<std::vector<ElemOp>( const std::vector<ElemOp> &in_f)> f) const
{
    if (!Terms.size()) return *this;
    Operator out_extra;
    Operator out_orig;
    Operator out(*this);
    for ( std::list<OpTerm>::iterator term_it = out.Terms.begin(); term_it != out.Terms.end(); term_it++) {
        OpTerm cur_term = *term_it;
        std::vector<ElemOp>& in_ops = cur_term.get<1>();  
        std::vector<ElemOp> out_ops = f(in_ops);
        //DEBUG(in_ops << "->" << out_ops);
        if (in_ops.size() != out_ops.size()) throw (exWrongOpSequence()); 
        unsigned int N=in_ops.size();
        // Here a hash check is done to check, that the result can be obtained by rearranging the operators
        std::size_t hash_in=0, hash_out=0;
        for (unsigned int i=0; i<N; ++i) {
            hash_in+=ElemOphash_value(in_ops[i]);
            hash_out+=ElemOphash_value(out_ops[i]);
            }
        assert (hash_in == hash_out );
        // Now time for moving the terms 
        for (unsigned int index_out = 0; index_out < out_ops.size(); index_out++) { // out_it loops over elements in resulting sequence in out_ops
            // Finds an iterator to the element corresponding to the one in out_it.
            std::vector<ElemOp>::iterator in_it = std::find(in_ops.begin()+index_out, in_ops.end(), out_ops[index_out]);
            unsigned int index_in = std::distance(in_ops.begin()+index_out, in_it)+index_out; 
            //DEBUG(index_out << " " << *in_it << " " << index_in);
            if (index_in - index_out > 0) 
                for (int index2 = int(index_in)-1; index2>=int(index_out); index2--) { // Here index2 can become negative, so comparison with an unsigned int can be passed.
                    std::pair<OpTerm,Operator> out_pair = elementary_swap_adjacent(cur_term, index2);
                    cur_term = out_pair.first;
                    out_extra+=out_pair.second;
                    }
            }
        out_orig+=cur_term;
    }
    out_extra = out_extra.rearrange(f);
    out_orig+=out_extra;
    out_orig.reduce();
    out_orig.prune();
    out_orig.sortTerms();
    return out_orig;
}

Operator Operator::getNormalOrdered() const
{
    return rearrange(NORMAL_ORDER);
}

bool Operator::operator==(const Operator &rhs)
{
    return ((this->getNormalOrdered().Terms) == (rhs.getNormalOrdered().Terms));
}

bool Operator::commutes(const Operator &rhs) const
{
    return ( (*this)*rhs == rhs*(*this));
}
boost::tuple<FockState,MelemType> Operator::actRight(const OpTerm &in, const FockState &ket)
{
    DEBUG(in << "|" << ket << ">");
    //ParticleIndex prev_pos_ = ket.size(); // Here we'll store the index of the last operator to speed up sign counting
    ParticleIndex prev_pos_ = 0;
    DEBUG(prev_pos_ <<"|" << ket);
    int sign=1;
    FockState bra = ket;
    MelemType Value; std::vector<ElemOp> in_ops;
    boost::tie(Value,in_ops)=in;
    unsigned int N=in_ops.size();
    for (int i=N-1; i>=0; i--) // Is the number of operator in OperatorSequence. Now we need to count them from back.
        {
            bool op; ParticleIndex ind;
            boost::tie(op,ind)=in_ops[i];
            DEBUG(op << "_" << ind);
            if (op == bra[ind] ) return boost::make_tuple(ERROR_FOCK_STATE, 0); // This is Pauli principle.
            if (ind > prev_pos_) 
                for (ParticleIndex j=prev_pos_; j<ind; ++j) { if (bra[j]) sign*=-1; } 
            else
                for (ParticleIndex j=prev_pos_; j>ind; j--) { if (bra[j]) sign*=-1; }
            bra[ind] = op; // This is c or c^+ acting
            //prev_pos_ = 0;
            
        }
    return boost::make_tuple(bra, Value*MelemType(sign));
}

template <class R>
inline bool __is_zero(std::pair<FockState,R> in){return (std::abs(in.second)<std::numeric_limits<RealType>::epsilon());};

std::map<FockState, MelemType> Operator::actRight(const FockState &ket) const
{
    DEBUG("!");
    std::map<FockState, MelemType> result1;
    for (std::list<OpTerm>::const_iterator it = Terms.begin(); it!=Terms.end(); it++)
        {
            FockState bra; 
            MelemType melem;
            boost::tie(bra,melem) = actRight(*it,ket);
            if (std::abs(melem)>1e-8) DEBUG(bra << "|*" << melem);
            if (bra!=ERROR_FOCK_STATE && std::abs(melem)>std::numeric_limits<RealType>::epsilon()) 
                result1[bra]+=melem;
        };
    // C++11 remove_if has a different behaviour, so this is a hck around. */
    std::map<FockState, MelemType> result2;
    std::remove_copy_if(result1.begin(), result1.end(), std::inserter(result2, result2.end()), __is_zero<MelemType> );//__is_zero<MelemType>);
    return result2;
}


MelemType Operator::getMatrixElement( const FockState & bra, const FockState &ket) const
{
    std::map<FockState, MelemType> output = this->actRight(ket);
    if (output.find(bra)==output.end()) 
        return 0;
    else { 
        return output[bra];
        }
}

MelemType Operator::getMatrixElement( const VectorType & bra, const VectorType &ket, const std::vector<FockState> &states) const
{
    if (bra.size()!=ket.size() || bra.size()!=states.size()) throw (exMelemVanishes()); 
    MelemType melem = 0.0;
    for (int i=0; i<ket.size(); ++i) {
        FockState current_state = states[i];
        MelemType overlap = ket[i];
        if (std::abs(overlap)>std::numeric_limits<RealType>::epsilon()) { 
            DEBUG(overlap << "," << current_state);
            std::map<FockState, MelemType> map1 = this->actRight(current_state);
            for (std::map<FockState, MelemType>::const_iterator it = map1.begin(); it!= map1.end(); it++) {
                FockState result_state = it->first;
                MelemType melem2 = it->second;
                DEBUG("\t<" << result_state << "|" << melem2);
                std::vector<FockState>::const_iterator it1 = std::find(states.begin(), states.end(), result_state);
                MelemType overlap2;
                if (it1 != states.end() ) { 
                    size_t j = std::distance(states.begin(), it1);
                #ifdef POMEROL_COMPLEX_MATRIX_ELEMENS
                    overlap2 = std::conj(bra(j));
                #else
                    overlap2 = bra(j);
                #endif
                    }
                else overlap2 = 0.0;
                //DEBUG(overlap2);
                //DEBUG("<" << result_state << "|" << overlap2 << "*" << melem << "*" << overlap << "|" << current_state << ">");
                melem += overlap2 * melem2 * overlap; 
                }   
            };  
        }; 
    return melem;
}

void Operator::reduce()
{
    int it1_pos=0;
    for (std::list<OpTerm>::iterator it1 = Terms.begin(); it1!=Terms.end(); it1++) {
        for (std::list<OpTerm>::iterator it2 = boost::next(it1); it2!=Terms.end();) {
            if (it2!=Terms.end()) {
                if (it2->get<1>() == it1->get<1>()) {
                    it1->get<0>()+=it2->get<0>();
                    it2 = Terms.erase(it2);
                    it1 = Terms.begin();
                    std::advance(it1,it1_pos);
                    }
                else it2++;
                }
            }
        it1_pos++;
    }
}

void Operator::prune(const RealType &Precision)
{
    for (std::list<OpTerm>::iterator it1 = Terms.begin(); it1!=Terms.end();) {
        if (std::abs(it1->get<0>()) < Precision) it1=Terms.erase(it1);
        else it1++;
        }
}

void Operator::sortTerms()
{
#ifdef POMEROL_COMPLEX_MATRIX_ELEMENS
    Terms.sort(__compareOpTerms);
#else
    Terms.sort();
#endif
}

const char* Operator::exWrongLabel::what() const throw(){
    return "Wrong labels";
};

const char* Operator::exWrongOpSequence::what() const throw(){
    std::stringstream s;
    s << "The term has wrong operator sequence!";
    return s.str().c_str();
};

const char* Operator::exMelemVanishes::what() const throw(){
    return "Matrix element vanishes"; 
};



Operator Operator::getCommutator(const Operator &rhs) const
{
    return (*this)*rhs - rhs*(*this);
}

Operator Operator::getAntiCommutator(const Operator &rhs) const
{
    return (*this)*rhs + rhs*(*this);
}


} // end of namespace Pomerol
