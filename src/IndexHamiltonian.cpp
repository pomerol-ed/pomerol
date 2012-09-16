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

/** \file hMatrix.cpp
**  \brief Implementation of IndexHamiltonian class - the Hamiltonian, writted in ParticleIndex Space.
** 
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "IndexHamiltonian.h"
#include <algorithm>
#include <sstream>

namespace Pomerol {

//
//IndexHamiltonian
//

IndexHamiltonian::IndexHamiltonian(const Lattice *L, const IndexClassification &IndexInfo):Operator(),L(L), IndexInfo(IndexInfo)
{
}

void IndexHamiltonian::prepare()
{
    // Read terms.
    for (unsigned int N=L->getTermStorage().getMaxTermOrder(); N; --N ) {
        if ( L->getTermStorage().getTerms(N).size())
        for (Lattice::TermList::const_iterator current=L->getTermStorage().getTerms(N).begin(); current!=L->getTermStorage().getTerms(N).end(); ++current) {
            std::vector<ParticleIndex> ind;
            ind.resize(N);
            for (unsigned int i=0; i<N; ++i) { 
                ind[i]=IndexInfo.getIndex((**current).SiteLabels[i], (**current).Orbitals[i], (**current).Spins[i]);
                };
            // Create a term out of the term in the lattice
            Operator::Term *T1 = new Operator::Term(N,(**current).OperatorSequence, ind, (**current).Value);
            mapTerms[N].push_back(T1);
            Terms->push_back(T1);
            } // end of Term loop
        } // end of for N
    //DEBUG("----------------------"); 
    // We now need to rearrange the terms to a given order for the easy access afterwards.
    for (unsigned int N=L->getTermStorage().getMaxTermOrder(); N; --N ) {
        std::map <unsigned int, std::list<Operator::Term*> >::iterator map_iterator = mapTerms.find(N);
        if (map_iterator!=mapTerms.end()) 
            for (std::list<Operator::Term*>::iterator it1=(map_iterator->second).begin(); it1!=(map_iterator->second).end(); it1++) {
                // get current term
                Operator::Term *T1 = *it1;
                // Rearrange it
                try {
                    boost::shared_ptr<std::list<Operator::Term*> > out=T1->rearrange(TERM_DEFAULT_SEQUENCE(N));
                    for (std::list<Operator::Term*>::iterator additional_terms = out->begin(); additional_terms != out->end(); additional_terms++) {
                        mapTerms[(**additional_terms).getN()].push_back(*additional_terms);
                        Terms->push_back(*additional_terms);
                        } // end of list iteration
                    }
                catch (Operator::Term::exWrongOpSequence)
                    {
                       INFO("Term " << *T1 << " couldn't be rearranged. ");// << TERM_DEFAULT_SEQUENCE(N) << " sequence.");
                    }
                } // end of map_iterator loop
            }; // end of terms order loop
};

const std::list<Operator::Term*> IndexHamiltonian::getTermsByOrder(unsigned int N) const
{
    return mapTerms.find(N)->second;
};

void IndexHamiltonian::printTerms(unsigned int N) const
{
    std::list<Operator::Term*> Temp = (this)->getTermsByOrder(N);
    for (std::list<Operator::Term*>::const_iterator it1=Temp.begin(); it1!=Temp.end(); ++it1) {
    INFO(**it1 );
    };

};

} // end of namespace Pomerol
