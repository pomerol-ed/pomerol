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

bool Operator::checkTerm(const OpTerm &in)
{
    if (std::abs(in.get<0>()) < std::numeric_limits<RealType>::epsilon()) return false;
    const std::vector<ElemOp> &ops = in.get<1>(); 
    unsigned int N=ops.size();
    for (unsigned int i=0; i<N; ++i) {  
        bool c; ParticleIndex ind;
        boost::tie(c,ind) = ops[i]; 
        int count_index=2*c-1; // This determines how many times current index Indices[i] is found. c^+ gives +1, c gives -1.
        //DEBUG(i << " " << Indices[i] << " " << count_index);
        for (unsigned int j=i+1; j<N; ++j) {
            bool c2; ParticleIndex ind2;
            boost::tie(c2,ind2)=ops[j];
            if (ind2 == ind ) {  //same indices
                count_index+=2*c2-1; 
                //DEBUG(j << " " << Indices[j] << " " << count_index);
                if ( count_index > 1 || count_index < -1 ) { ERROR("This term vanishes. "); return false; }; 
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
    Terms.reset( new std::list<OpTerm> );
}

Operator::Operator(boost::shared_ptr<std::list<OpTerm> > Terms) : Terms(Terms)
{
    Terms->remove_if(this->checkTerm);
}

Operator::Operator(const OpTerm& term)
{
    Terms.reset( new std::list<OpTerm> );
    if (checkTerm(term)) Terms->push_back(term);
}

boost::shared_ptr<std::list<OpTerm> > Operator::getTerms() const
{
    return Terms;
}

Operator::~Operator()
{
    Terms.reset();
}
bool Operator::isEmpty() const
{
    return (Terms->size()==0);
}

std::ostream& operator<< (std::ostream& output, const Operator& out)
{
    for (std::list<OpTerm>::const_iterator it = out.Terms->begin(); it!=out.Terms->end(); it++) {
        if (it!=out.Terms->begin())  output << " + ";
        output << *it;
        //output << it->get<0>();
        };
    return output;
}

void Operator::printAllTerms() const
{
    INFO(*this);
}

Operator& Operator::operator+=(const Operator &rhs)
{
    if ( rhs.Terms->size() ) std::copy_backward(rhs.Terms->begin(), rhs.Terms->end(), Terms->end()); 
    return *this;
}

void Operator::add(const OpTerm &rhs)
{
    if (checkTerm(rhs)) Terms->push_back(rhs);
}

Operator& Operator::operator+=(const OpTerm &rhs)
{
    if (checkTerm(rhs)) Terms->push_back(rhs);
    return *this;
}

const Operator Operator::operator+(const Operator &rhs) const
{
    Operator out(*this);
    out+=rhs;
    return out;
}

Operator Operator::elementary_swap_adjacent(OpTerm &in, unsigned int position, bool force_ignore_commutation)
{
    Operator out;
    std::vector<ElemOp>& in_ops = in.get<1>();  MelemType Value=in.get<0>();
    if ( in_ops[position].get<1>() != in_ops[position+1].get<1>() || force_ignore_commutation ) {
        in.get<0>()*=(-1.); 
        boost::swap(in_ops[position], in_ops[position+1]);
        }
    else {
        if (in_ops[position].get<0>() == in_ops[position+1].get<0>()) throw (exWrongOpSequence()); 
        OpTerm term;
        std::copy_backward(in_ops.begin(), in_ops.begin()+position, term.get<1>().end());
        std::copy_backward(in_ops.begin()+position+1, in_ops.end(), term.get<1>().end());
        term.get<0>() = Value;
        out+=term;
        elementary_swap(in, position, true);
         };
    return out;
}

Operator Operator::elementary_swap(OpTerm &in, unsigned int position1, unsigned int position2, bool force_ignore_commutation)
{
    Operator out;
    std::vector<ElemOp> in_ops;  MelemType Value;
    boost::tie(Value,in_ops)=in;
    if (position2 == position1) return out;
    if (position2 < position1) std::swap(position2, position1);
    for (unsigned int i=position1; i<position2; i++) out+=elementary_swap_adjacent(in, i);
    for (unsigned int i=position2-2; i>=position1; i++) out+=elementary_swap_adjacent(in, i);
    return out;
}

void Operator::rearrange(boost::function<std::vector<ElemOp>( const std::vector<ElemOp> &in_f)> f)
{
    if (!Terms->size()) return;
    Operator out;
    for ( std::list<OpTerm>::iterator term_it = Terms->begin(); term_it != Terms->end(); term_it++) {
        std::vector<ElemOp>& in_ops = term_it->get<1>();  
        std::vector<ElemOp> out_ops = f(in_ops);
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
            unsigned int index_in = std::distance(in_ops.begin(), in_it); 
            if (index_in - index_out > 0) for (unsigned int index2 = index_in-1; index2>=index_out; index2--) out+=elementary_swap_adjacent(*term_it,index2);
        }
    }
    out.rearrange(f);
    (*this)+=out;
    this->reduce();
    this->prune();
}

void Operator::makeNormalOrder()
{
    return rearrange(NORMAL_ORDER);
}

boost::tuple<FockState,MelemType> Operator::actRight(const OpTerm &in, const FockState &ket)
{
    ParticleIndex prev_pos_ = 0; // Here we'll store the index of the last operator to speed up sign counting
    int sign=1;
    FockState bra = ket;
    MelemType Value; std::vector<ElemOp> in_ops;
    boost::tie(Value,in_ops)=in;
    unsigned int N=in_ops.size();
    for (int i=N-1; i>=0; i--) // Is the number of operator in OperatorSequence. Now we need to count them from back.
        {
            bool op; ParticleIndex ind;
            boost::tie(op,ind)=in_ops[i];
            if (op == bra[ind] ) return boost::make_tuple(ERROR_FOCK_STATE, 0); // This is Pauli principle.
            bra[ind] = op; // This is c or c^+ acting
            if (ind > prev_pos_) 
                for (ParticleIndex j=prev_pos_; j<ind; ++j) { if (ket[j]) sign*=-1; } 
            else
                for (ParticleIndex j=prev_pos_; j>ind; j--) { if (ket[j]) sign*=-1; }
            
        }
    return boost::make_tuple(bra, Value*MelemType(sign));
}

std::map<FockState, MelemType> Operator::actRight(const FockState &ket) const
{
    std::map<FockState, MelemType> result1;
    for (std::list<OpTerm>::const_iterator it = Terms->begin(); it!=Terms->end(); it++)
        {
            FockState bra; 
            MelemType melem;
            boost::tie(bra,melem) = actRight(*it,ket);
            if (bra!=ERROR_FOCK_STATE && std::abs(melem)>std::numeric_limits<RealType>::epsilon()) 
                result1[bra]+=melem;
        }
    for (std::map<FockState, MelemType>::iterator it1 = result1.begin(); it1!=result1.end(); it1++) 
        if ( std::abs(it1->second)<std::numeric_limits<RealType>::epsilon() ) result1.erase(it1); 
    return result1;
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


void Operator::reduce()
{
    int it1_pos=0;
    for (std::list<OpTerm>::iterator it1 = Terms->begin(); it1!=Terms->end(); it1++) {
        for (std::list<OpTerm>::iterator it2 = boost::next(it1); it2!=Terms->end();) {
            if (it2!=Terms->end()) {
                if (it2->get<1>() == it1->get<1>()) {
                    it1->get<0>()+=it2->get<0>();
                    it2 = Terms->erase(it2);
                    it1 = Terms->begin();
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
    for (std::list<OpTerm>::iterator it1 = Terms->begin(); it1!=Terms->end(); it1++)
        if (std::abs(it1->get<0>()) < Precision) it1=Terms->erase(it1);
}

const char* Operator::exWrongLabel::what() const throw(){
    return "Wrong labels";
};

const char* Operator::exWrongOpSequence::what() const throw(){
    std::stringstream s;
    s << "The term has wrong operator sequence!";
    return s.str().c_str();
};


/*
void Operator::makeNormalOrder()
{
    std::list<OpTerm*> output; // Here we will store the normal ordered output terms.
    for (std::list<OpTerm*>::const_iterator it = Terms->begin(); it!=Terms->end(); it++) {
        boost::shared_ptr<std::list<OpTerm* > > out = (**it).makeNormalOrder();
        output.splice(boost::prior(output.end()),*out);
    }
    Terms->splice(boost::prior(Terms->end()),output);
    this->reduce();
    this->prune();
}



*/


/*
Operator Operator::getCommutator(const Operator &rhs) const
{
   // DEBUG("!!" << rhs.Terms->size());
    boost::shared_ptr<std::list<OpTerm*> > output ( new std::list<OpTerm*>);
    for (std::list<OpTerm*>::const_iterator it = Terms->begin(); it!=Terms->end(); it++)
        for (std::list<OpTerm*>::const_iterator it2 = rhs.Terms->begin(); it2!=rhs.Terms->end(); it2++) {
   //         DEBUG("Commuting" << **it << " and " << **it2);
            boost::shared_ptr<std::list<OpTerm* > > out2 = (**it).getCommutator(**it2); 
    //        DEBUG(out2->size() << " terms generated.");
            output->splice(boost::prior(output->end()),*out2);
    //        DEBUG(output->size() << " in output.");
    }
    Operator out(output);
    return out;
}

bool Operator::commutes(const Operator &rhs) const
{
    bool out;
    Operator out_return(this->getCommutator(rhs));
    INFO("+++++++");
    out_return.printAllTerms();
    INFO("+++++++");
    out_return.makeNormalOrder();
    out_return.printAllTerms();
    INFO("+++++++");
    out_return.reduce();
    out_return.prune();
    out_return.printAllTerms();
//    DEBUG(out_return.Terms->size());
    return out;
}

*/

//
//Operator::Term
//

/*
boost::shared_ptr<std::list<Operator::Term*> > Operator::Term::rearrange(const std::vector<bool> & DesiredSequence)
{
    if (DesiredSequence.size() != OperatorSequence.size() ) throw (exWrongOpSequence());
    boost::shared_ptr<std::list<Operator::Term*> > out ( new std::list<Operator::Term*> );
    if (OperatorSequence == DesiredSequence) { return out; } // Nothing is needed to do then.

    for (unsigned int i=0; i<N-1; ++i) { 
        if (OperatorSequence[i]!=DesiredSequence[i]) {
            unsigned int j=i+1;
            for (; j<N && ( OperatorSequence[j] == OperatorSequence[i] || OperatorSequence[j] == DesiredSequence[j]); ++j) {};  // finding element to change
            if (j==N) throw (exWrongOpSequence()); // exit if there is no way of doing the rearrangement.
            if (N==2) { elementary_swap(0,true); return out; }; // If there are only two operators - just swap them and go away.

            // Now check if any operations with anticommutation relation is required.
            bool needNewTerms=false;
            for (unsigned int k=i+1; k<j && !needNewTerms; k++) needNewTerms = (Indices[k]==Indices[i] || Indices[k]==Indices[j]);

            if ( !needNewTerms ) { // Operators anticommute then. Constant energy levels are omitted.
                    Value*=(-1.); 
                    OperatorSequence[i]=!OperatorSequence[i]; 
                    OperatorSequence[j]=!OperatorSequence[j]; 
                    std::swap(Indices[i], Indices[j] );
                    }
            else { // Operators do not anticommute - swap will construct additional term.
                for (unsigned int k=j-1; k>=i; k--) { // move an operator at position j to the left to the position i.
                    std::list<Operator::Term*> out_temp = *(elementary_swap(k)); 
                    out->resize(out->size()+out_temp.size());
                    std::copy_backward(out_temp.begin(), out_temp.end(), out->end()); 
                    if (k==0) break; // exit, since unsigned int loop is done.
                    };
                for (unsigned int k=i+1; k<j; k++) { // move the operator at position i+1 to the right to the position j.
                    std::list<Operator::Term*> out_temp = *(elementary_swap(k)); 
                    out->resize(out->size()+out_temp.size());
                    std::copy_backward(out_temp.begin(), out_temp.end(), out->end()); 
                    };
                }; // end of else
            }; // end of element check
        } // end of i loop
    // Now we have the right order of booleans - what about indices? Indices should be rearranged prior to this
    return out;
}

void Operator::Term::reorder(bool ascend)
{
    //DEBUG("reordering " << *this);
    if (!N) return;
    assert ( OperatorSequence[0]==1 && OperatorSequence[N-1]==0);
    for (unsigned int i=0; i<N/2; ++i)
        for (unsigned int j=i*(!ascend); j<N/2-i*ascend-1; ++j)
        {
            if (ascend) {
                if (Indices[j+1] < Indices[j]) elementary_swap(j,true);
                if (Indices[j+1+N/2] < Indices[j+N/2]) elementary_swap(j+N/2,true);
                }
            else {
                if (Indices[j+1] > Indices[j]) elementary_swap(j,true);
                if (Indices[j+1+N/2] > Indices[j+N/2]) elementary_swap(j+N/2,true); // indices are unsigned int, so this is stable
            }
        }
}

boost::shared_ptr<std::list<Operator::Term*> > Operator::Term::makeNormalOrder()
{
    //boost::shared_ptr<std::list<Operator::Term*> > out ( new std::list<Operator::Term*> );
    //return out;
    std::vector<bool> normalOrderedSequence;
    for (ParticleIndex i=0; i<N; ++i) {
        if ( !OperatorSequence[i] ) normalOrderedSequence.push_back(0);
        else normalOrderedSequence.insert(normalOrderedSequence.begin(),1);
        //else normalOrderedSequence.resize(normalOrderedSequence.size()+1,1); // for a bitset
        }
    boost::shared_ptr<std::list<Operator::Term*> > out = this->rearrange(normalOrderedSequence);
    //DEBUG("Main: " << *this);
    this->reorder();
    for (std::list<Operator::Term*>::iterator iter_it = out->begin(); iter_it!=out->end(); ++iter_it) {
            //DEBUG("Additional: " << **iter_it);
            boost::shared_ptr<std::list<Operator::Term*> > out2 = (*iter_it)->makeNormalOrder();
            out->splice(boost::prior(out->end()),*out2);
        }
    return out;
}

MelemType Operator::Term::getMatrixElement( const FockState & bra, const FockState &ket){
    MelemType result;
    FockState bra2;
    boost::tie(bra2, result) = this->actRight(ket);
    return (bra2 == bra)?result:0;
}


unsigned int Operator::Term::getN(){
    return N;
}

std::ostream& operator<< (std::ostream& output, const Operator::Term& out)
{
    output << out.Value << "*"; 
    for (unsigned int i=0; i<out.N; ++i) output << ((out.OperatorSequence[i])?"c^{+}":"c") << "_" << out.Indices[i];
    return output; 
}

bool Operator::Term::isExactlyEqual(const Operator::Term &rhs) const
{
    bool out=(N==rhs.N && Value==rhs.Value);
    if (!out) return false;
    for (unsigned int i=0; i<N; ++i) { out = (out && OperatorSequence[i] == rhs.OperatorSequence[i] && Indices[i] == rhs.Indices[i]); };
    return out;
}

bool Operator::Term::operator==(const Operator::Term &rhs) const
{
    bool out=(rhs.isExactlyEqual(*this));
    if (out) return true;
    Operator::Term this_copy(*this);
    Operator::Term rhs_copy(rhs);
    boost::shared_ptr<std::list<Operator::Term*> > list_lhs = this_copy.makeNormalOrder();
    boost::shared_ptr<std::list<Operator::Term*> > list_rhs = rhs_copy.makeNormalOrder();
    reduce(list_lhs);
    reduce(list_rhs);
//    DEBUG(this_copy << "," << rhs_copy);
//    DEBUG(list_lhs->size());
    if (!(this_copy.isExactlyEqual(rhs_copy)) || list_lhs->size() != list_rhs->size() ) return false;
    out = true;
    for ( std::pair<std::list<Operator::Term*>::iterator, std::list<Operator::Term*>::iterator> pair_it = std::make_pair(list_lhs->begin(), list_rhs->begin()); 
          pair_it.first!=list_lhs->end() && pair_it.second!=list_rhs->end(); ) { 
//        DEBUG(**pair_it.first);
//        DEBUG(**pair_it.second);
        out = (out && (**pair_it.first).isExactlyEqual(**pair_it.second));
        pair_it.first++; pair_it.second++;
        }
    return out;
}



boost::shared_ptr<std::list<Operator::Term*> > Operator::Term::getCommutator(const Operator::Term &rhs) const
{
    boost::shared_ptr<std::list<Operator::Term*> > out ( new std::list<Operator::Term*> );
    int Ntotal = N+rhs.N;
    
    std::vector<bool> Seq2(Ntotal);
    std::vector<unsigned int> Ind2(Ntotal);

    for (unsigned int i=0; i<N; ++i)      { Seq2[i] = OperatorSequence[i]; Ind2[i] = Indices[i]; };
    for (unsigned int i=N; i<Ntotal; ++i) { Seq2[i] = rhs.OperatorSequence[i-N]; Ind2[i] = rhs.Indices[i-N]; };
    try {
            out->push_back(new Term(Ntotal, Seq2, Ind2, Value));
        }
    catch (Operator::Term::exWrongOpSequence &e)
        {
            ERROR("Commutator [" << *this << "," << rhs << "] creates a vanishing term");
        }

    for (unsigned int i=0; i<rhs.N; ++i)      { Seq2[i] = rhs.OperatorSequence[i]; Ind2[i] = rhs.Indices[i]; };
    for (unsigned int i=rhs.N; i<Ntotal; ++i) { Seq2[i] = OperatorSequence[i-rhs.N]; Ind2[i] = Indices[i-rhs.N]; };
    try {
            out->push_back(new Term(Ntotal, Seq2, Ind2, -Value));
        }
    catch (Operator::Term::exWrongOpSequence &e)
        {
            ERROR("Commutator [" << *this << "," << rhs << "] creates a vanishing term");
        };

    return out;
}

bool Operator::Term::commutes(const Operator::Term &rhs) const
{
    boost::shared_ptr<std::list<Operator::Term*> > out_terms = this->getCommutator(rhs);
    if (out_terms->size() == 0) return true;
    if (out_terms->size() == 1) return false;
    assert (out_terms->size() == 2);
    Operator::Term* term1 = (*out_terms->begin());
    Operator::Term* term2 = (*out_terms->begin()++);
    return (*term1 == *term2);
}

void Operator::Term::prune(boost::shared_ptr<std::list<Operator::Term*> > Terms, const RealType &Precision)
{
    //std::remove_if (Terms->begin(), Terms->end(), std::abs(boost::bind(boost::lambda::_1)->Value) < Precision);
    for (std::list<Operator::Term*>::iterator it1 = Terms->begin(); it1!=Terms->end(); it1++)
        if (std::abs((**it1).Value) < Precision) it1=Terms->erase(it1);
}
*/

} // end of namespace Pomerol
