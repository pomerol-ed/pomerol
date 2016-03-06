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

#include "pomerol/IndexHamiltonian.h"
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
            Operator tmp; 
            for (unsigned int i=0; i<N; ++i) { 
                ParticleIndex i1 = IndexInfo.getIndex((**current).SiteLabels[i], (**current).Orbitals[i], (**current).Spins[i]);
                Operator t1 = ((**current).OperatorSequence[i]==Lattice::Term::creation)?OperatorPresets::c_dag(i1):OperatorPresets::c(i1);
                if (tmp.isEmpty()) tmp=t1;
                else tmp*=t1; 
                };
            // Create a term out of the term in the lattice
            
            (*this)+=(**current).Value*tmp;
            } // end of Term loop
        } // end of for N
};

/*
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
*/

} // end of namespace Pomerol
