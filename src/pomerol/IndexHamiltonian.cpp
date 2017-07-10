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
