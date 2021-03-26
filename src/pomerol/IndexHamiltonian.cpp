#include "pomerol/IndexHamiltonian.h"
#include <algorithm>
#include <sstream>

namespace Pomerol {

//
//IndexHamiltonian
//

template<bool Complex>
IndexHamiltonian<Complex>::IndexHamiltonian(const Lattice<Complex> *L,
                                            const IndexClassification<Complex> &IndexInfo):
  Operator<Complex>(), L(L), IndexInfo(IndexInfo)
{
}

template<bool Complex>
void IndexHamiltonian<Complex>::prepare()
{
    // Read terms.
    for (unsigned int N=L->getTermStorage().getMaxTermOrder(); N; --N ) {
        if ( L->getTermStorage().getTerms(N).size())
        for (auto current=L->getTermStorage().getTerms(N).begin(); current!=L->getTermStorage().getTerms(N).end(); ++current) {
            Operator<Complex> tmp;
            for (unsigned int i=0; i<N; ++i) {
                ParticleIndex i1 = IndexInfo.getIndex((**current).SiteLabels[i], (**current).Orbitals[i], (**current).Spins[i]);
                Operator<Complex> t1 = ((**current).OperatorSequence[i]==Lattice<Complex>::Term::creation)
                                       ? OperatorPresets::c_dag<Complex>(i1) : OperatorPresets::c<Complex>(i1);
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

template class IndexHamiltonian<false>;
template class IndexHamiltonian<true>;

} // end of namespace Pomerol
