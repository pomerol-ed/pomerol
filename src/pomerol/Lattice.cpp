#include "pomerol/Lattice.h"
#include "pomerol/LatticePresets.h"
#include <fstream>
#include <algorithm>

namespace Pomerol{

//
// Site
//

Site::Site(const std::string& Label, unsigned short OrbitalSize, unsigned short SpinSize):Label(Label), OrbitalSize(OrbitalSize), SpinSize(SpinSize)
{
};


std::ostream& operator<<(std::ostream& output, const Site& out)
{
    output << "Site \"" << out.Label << "\", " << out.OrbitalSize << " orbital" << ((out.OrbitalSize>1)?"s":"") << ", " << out.SpinSize << " spin" << ((out.SpinSize>1)?"s":"") << ".";
	return output;
}

//
// Lattice::Term
//

template<bool Complex>
Lattice<Complex>::Term::Term (unsigned int N):N(N)
{
    OperatorSequence.resize(N);
    SiteLabels.resize(N);
    Spins.resize(N);
    Orbitals.resize(N);
    for (unsigned int i=0; i<N; ++i) { OperatorSequence[i]=false; SiteLabels[i]=""; Spins[i]=0; Orbitals[i]=0.0; };
    Value=0.0;
};

template<bool Complex>
Lattice<Complex>::Term::Term(unsigned int N,
                             bool * OperatorSequence_,
                             MelemType<Complex> Value_,
                             std::string * SiteLabels_,
                             unsigned short * Orbitals_,
                             unsigned short *Spins_):
N(N)
{
  OperatorSequence.assign( OperatorSequence_, OperatorSequence_+N );
  SiteLabels.assign( SiteLabels_, SiteLabels_+N );
  Spins.assign( Spins_, Spins_+N );
  Orbitals.assign( Orbitals_, Orbitals_+N );
  Value=Value_;
}

template<bool Complex>
Lattice<Complex>::Term::Term (const Lattice::Term &in):N(in.N), OperatorSequence(in.OperatorSequence), SiteLabels(in.SiteLabels), Spins(in.Spins), Orbitals(in.Orbitals), Value(in.Value)
{
};

template<bool Complex>
unsigned int Lattice<Complex>::Term::getOrder() const { return N; };

template<bool Complex>
std::ostream& operator<< (std::ostream& output, const typename Lattice<Complex>::Term& out)
{
    output << out.Value << "*";
    for (unsigned int i=0; i<out.N; ++i) output << ((out.OperatorSequence[i])?"c^{+}":"c") << "_{" << out.SiteLabels[i] << "," << out.Orbitals[i] << "," << out.Spins[i] << "}" ;
    return output;
};

//
// Lattice::TermStorage
//

template<bool Complex>
Lattice<Complex>::TermStorage::TermStorage()
{
    MaxTermOrder=0;
};

template<bool Complex>
int Lattice<Complex>::TermStorage::addTerm(const Lattice::Term *T)
{
    unsigned int N = T->getOrder();
    Terms[N].push_back(new Term(*T));
    MaxTermOrder=(MaxTermOrder<N)?N:MaxTermOrder;
    return 0;
};

template<bool Complex>
const unsigned int Lattice<Complex>::TermStorage::getMaxTermOrder() const
{
    return MaxTermOrder;
}

template<bool Complex>
const typename Lattice<Complex>::TermList &Lattice<Complex>::TermStorage::getTerms (unsigned int N) const
{
   auto it1 = Terms.find(N);
    if (Terms.find(N)!=Terms.end())
        {
            return it1->second;
        }
    else return *(new TermList ());
};

//
// Lattice
//

template<bool Complex>
Lattice<Complex>::Lattice():Terms(new TermStorage)
{
};

template<bool Complex>
Lattice<Complex>::~Lattice(){
delete Terms;
};

template<bool Complex>
Lattice<Complex>::Lattice(const Lattice &l) : Sites(l.Sites) {
 Terms = new TermStorage(*l.Terms);
}

template<bool Complex>
auto Lattice<Complex>::getSiteMap() const -> const SiteMap&
{
    return Sites;
}

template<bool Complex>
auto Lattice<Complex>::getTermStorage() const -> const TermStorage&
{
    return *Terms;
}

template<bool Complex>
void Lattice<Complex>::printTerms(unsigned int n)
{
TermList Temp = Terms->getTerms(n);
    for (auto it1 = Temp.begin(); it1!=Temp.end(); ++it1) {
        INFO(**it1 );
    };
}

template<bool Complex>
void Lattice<Complex>::printSites() const
{
    for (SiteMap::const_iterator it1=Sites.begin(); it1!=Sites.end(); ++it1) {
            INFO(*(it1->second));
        };
}

template<bool Complex>
void Lattice<Complex>::addSite(Site* S)
{
    Sites[S->Label]= S ;
}

template<bool Complex>
void Lattice<Complex>::addSite(const std::string &Label, unsigned short orbitals, unsigned short spins)
{
    addSite(new Site(Label, orbitals, spins));
}

template<bool Complex>
void Lattice<Complex>::addTerm(const Lattice::Term *T)
{
    unsigned int N=T->getOrder();
    for (unsigned int i=0; i<N; ++i) { // some checks to avoid bad terms.
        if (Sites.find(T->SiteLabels[i]) == Sites.end()) { ERROR("Sites, specified in this Term " << *T << " do not exist."); throw (exWrongLabel()); };
        if (T->Orbitals[i]>=Sites[T->SiteLabels[i]]->OrbitalSize) { ERROR("Wrong orbital indices in term " << *T << "."); throw (exWrongLabel()); };
        if (T->Spins[i]>=Sites[T->SiteLabels[i]]->SpinSize) { ERROR("Wrong spin indices in term " << *T << "."); throw (exWrongLabel()); };
    };
    if ( std::abs(T->Value) ) Terms->addTerm(T);
}

template<bool Complex>
const Site& Lattice<Complex>::getSite(const std::string& Label) const
{
    std::map<std::string, Site*>::const_iterator it1=Sites.find(Label);
    if (it1!=Sites.end()) throw (exWrongLabel());
    return *(it1->second);
}

template<bool Complex>
const char* Lattice<Complex>::exWrongLabel::what() const throw(){
    return "Wrong requested Label";
};

template class Lattice<false>;
template class Lattice<true>;

} // end of namespace Pomerol
