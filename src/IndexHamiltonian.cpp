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

/** \file hMatrix.cpp
**  \brief Implementation of IndexHamiltonian class - the Hamiltonian, writted in ParticleIndex Space.
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#include "IndexHamiltonian.h"
#include <algorithm>

namespace Pomerol {

//
//IndexHamiltonian::Term
//

IndexHamiltonian::Term::Term (const unsigned int Order, const std::vector<bool>&  Sequence, const std::vector<ParticleIndex> & Indices, RealType Value):
    Order(Order), OperatorSequence(Sequence), Indices(Indices), Value(Value)
{
    if (Sequence.size()!=Order || Indices.size()!=Order) throw(exWrongLabel());
    for (unsigned int i=0; i<Order; ++i) 
        for (unsigned int j=i+1; j<Order; ++j)
            if (Indices[i]==Indices[j] && Sequence[i] == Sequence[j] ) { DEBUG(Sequence[i] << "," << Sequence[j]); ERROR("Cannot have two same operators in a term. "); throw (exWrongOpSequence()); }; 
}

std::list<IndexHamiltonian::Term*> IndexHamiltonian::Term::rearrange(const std::vector<bool> & DesiredSequence)
{
    if (DesiredSequence.size() != OperatorSequence.size() ) throw (exWrongOpSequence());
    std::list<IndexHamiltonian::Term*> out;
    if (OperatorSequence == DesiredSequence) { return out; } // Nothing is needed to do then.

    for (unsigned int i=0; i<Order-1; ++i) { 
        if (OperatorSequence[i]!=DesiredSequence[i]) {
            unsigned int j=i+1;
            for (; j<Order && ( OperatorSequence[j] == OperatorSequence[i] || OperatorSequence[j] == DesiredSequence[j]); ++j) {};  // finding element to change
            if (j==Order) throw (exWrongOpSequence()); // exit if there is no way of doing the rearrangement.
            if (Order==2) { elementary_swap(0,true); return out; }; // If there are only two operators - just swap them and go away.

            // Now check if any operations with anticommutation relation is required.
            bool needNewTerms=false;
            for (unsigned int k=i+1; k<j && !needNewTerms; k++) needNewTerms = (Indices[k]==Indices[i] || Indices[k]==Indices[j]);

            if ( !needNewTerms ) { // Operators anticommute then. Constant energy levels are omitted.
                    Value*=(-1.); 
                    OperatorSequence[i]=!OperatorSequence[i]; 
                    OperatorSequence[j]=!OperatorSequence[j]; 
                    std::swap(Indices[i], Indices[j] );
                    }
            else { // Operators do not commute - swap will construct additional term.
                for (unsigned int k=j-1; k>=i; k--) { // move an operator at position j to the left to the position i.
                    std::list<IndexHamiltonian::Term*> out_temp = elementary_swap(k); 
                    out.resize(out.size()+out_temp.size());
                    std::copy_backward(out_temp.begin(), out_temp.end(), out.end()); 
                    if (k==0) break; // exit, since unsigned int loop is done.
                    };
                for (unsigned int k=i+1; k<j; k++) { // move the operator at position i+1 to the right to the position j.
                    std::list<IndexHamiltonian::Term*> out_temp = elementary_swap(k); 
                    out.resize(out.size()+out_temp.size());
                    std::copy_backward(out_temp.begin(), out_temp.end(), out.end()); 
                    };
                }; // end of else
            }; // end of element check
        } // end of i loop
    return out;
}

std::list<IndexHamiltonian::Term*> IndexHamiltonian::Term::elementary_swap(unsigned int position, bool force_ignore_commutation)
{
    std::list<IndexHamiltonian::Term*> out;
    if ( Indices[position] != Indices[position+1] || force_ignore_commutation ) {
        Value*=(-1.); 
        bool tmp = OperatorSequence[position];
        OperatorSequence[position]=OperatorSequence[position+1];
        OperatorSequence[position+1]=tmp;
        std::swap(Indices[position], Indices[position+1] );
        }
    else {
        std::vector<bool> Seq2(Order-2);
        std::vector<unsigned int> Ind2(Order-2);
        for (unsigned int i=0; i<position; ++i) { Seq2[i] = OperatorSequence[i]; Ind2[i] = Indices[i]; };
        for (unsigned int i=position+2; i<Order; ++i) { Seq2[i-2] = OperatorSequence[i]; Ind2[i-2] = Indices[i]; };
        out.push_back(new Term(Order-2, Seq2, Ind2, Value));
        elementary_swap(position, true);
         };
    return out;
}

const char* IndexHamiltonian::Term::exWrongLabel::what() const throw(){
    return "Wrong labels";
};

const char* IndexHamiltonian::Term::exWrongOpSequence::what() const throw(){
    return "The term can not be rearranged!";
};

std::ostream& operator<< (std::ostream& output, const IndexHamiltonian::Term& out)
{
    output << out.Value << "*"; 
    for (unsigned int i=0; i<out.Order; ++i) output << ((out.OperatorSequence[i])?"c^{+}":"c") << "_" << out.Indices[i];
    return output; 
}

//
//IndexHamiltonian
//

IndexHamiltonian::IndexHamiltonian(const Lattice *L, const IndexClassification &IndexInfo):L(L), IndexInfo(IndexInfo)
{
}

void IndexHamiltonian::prepare()
{
    // Read terms.
    DEBUG(L->getTermStorage().getMaxTermOrder());
    for (unsigned int N=L->getTermStorage().getMaxTermOrder(); N; --N ) {
        if ( L->getTermStorage().getTerms(N).size())
        for (Lattice::TermList::const_iterator current=L->getTermStorage().getTerms(N).begin(); current!=L->getTermStorage().getTerms(N).end(); ++current) {
            std::vector<ParticleIndex> ind;
            ind.resize(N);
            for (unsigned int i=0; i<N; ++i) { 
                ind[i]=IndexInfo.getIndex((**current).SiteLabels[i], (**current).Orbitals[i], (**current).Spins[i]);
                };
            // Create a term out of the term in the lattice
            IndexHamiltonian::Term *T1 = new IndexHamiltonian::Term(N,(**current).OperatorSequence, ind, (**current).Value);
            Terms[N].push_back(T1);
            } // end of Term loop
        } // end of for N
    //DEBUG("----------------------"); 
    // We now need to rearrange the terms to a given order for the easy access afterwards.
    for (unsigned int N=L->getTermStorage().getMaxTermOrder(); N; --N ) {
        std::map <unsigned int, std::list<Term*> >::iterator map_iterator = Terms.find(N);
        if (map_iterator!=Terms.end()) 
            for (std::list<IndexHamiltonian::Term*>::iterator it1=(map_iterator->second).begin(); it1!=(map_iterator->second).end(); it1++) {
                // get current term
                IndexHamiltonian::Term *T1 = *it1;
                // Rearrange it
                std::list<IndexHamiltonian::Term*> out=T1->rearrange(TERM_DEFAULT_SEQUENCE(N));
       //         DEBUG(*T1);
                for (std::list<IndexHamiltonian::Term*>::iterator additional_terms = out.begin(); additional_terms != out.end(); additional_terms++) {
                        Terms[(**additional_terms).Order].push_back(*additional_terms);
                    } // end of list iteration
                } // end of map_iterator loop
            }; // end of terms order loop
};

const std::list<IndexHamiltonian::Term*> IndexHamiltonian::getTerms(unsigned int N) const
{
    return Terms.find(N)->second;
};

} // end of namespace Pomerol
