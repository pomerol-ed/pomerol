#include "pomerol/Lattice.h"
#include "pomerol/LatticePresets.h"
#include <fstream>
#include <algorithm>

namespace Pomerol{

//
// Lattice::Site
//

Lattice::Site::Site(const std::string& Label, unsigned short OrbitalSize, unsigned short SpinSize):Label(Label), OrbitalSize(OrbitalSize), SpinSize(SpinSize)
{
};


std::ostream& operator<<(std::ostream& output, const Lattice::Site& out)
{
    output << "Site \"" << out.Label << "\", " << out.OrbitalSize << " orbital" << ((out.OrbitalSize>1)?"s":"") << ", " << out.SpinSize << " spin" << ((out.SpinSize>1)?"s":"") << ".";
	return output;
}

//
// Lattice::Term
//

Lattice::Term::Term (unsigned int N):N(N)
{
    OperatorSequence.resize(N);
    SiteLabels.resize(N);
    Spins.resize(N);
    Orbitals.resize(N);
    for (unsigned int i=0; i<N; ++i) { OperatorSequence[i]=false; SiteLabels[i]=""; Spins[i]=0; Orbitals[i]=0.0; }; 
    Value=0.0; 
};


Lattice::Term::Term(unsigned int N, bool * OperatorSequence_, MelemType Value_, std::string * SiteLabels_, unsigned short * Orbitals_, unsigned short *Spins_):
N(N)
{
  OperatorSequence.assign( OperatorSequence_, OperatorSequence_+N );
  SiteLabels.assign( SiteLabels_, SiteLabels_+N );
  Spins.assign( Spins_, Spins_+N );
  Orbitals.assign( Orbitals_, Orbitals_+N );
  Value=Value_;
}

Lattice::Term::Term (const Lattice::Term &in):N(in.N), OperatorSequence(in.OperatorSequence), SiteLabels(in.SiteLabels), Spins(in.Spins), Orbitals(in.Orbitals), Value(in.Value)
{
};
unsigned int Lattice::Term::getOrder() const { return N; };

std::ostream& operator<< (std::ostream& output, const Lattice::Term& out)
{   
    output << out.Value << "*"; 
    for (unsigned int i=0; i<out.N; ++i) output << ((out.OperatorSequence[i])?"c^{+}":"c") << "_{" << out.SiteLabels[i] << "," << out.Orbitals[i] << "," << out.Spins[i] << "}" ; 
    return output; 
};

//
// Lattice::TermStorage
//

Lattice::TermStorage::TermStorage()
{
    MaxTermOrder=0;
};

int Lattice::TermStorage::addTerm(const Lattice::Term *T)
{
    unsigned int N = T->getOrder();
    Terms[N].push_back(new Term(*T));
    MaxTermOrder=(MaxTermOrder<N)?N:MaxTermOrder;
    return 0;
};

const unsigned int Lattice::TermStorage::getMaxTermOrder() const
{
    return MaxTermOrder;  
}

const Lattice::TermList &Lattice::TermStorage::getTerms (unsigned int N) const
{
   std::map<unsigned int, Lattice::TermList>::const_iterator it1=Terms.find(N); 
    if (Terms.find(N)!=Terms.end()) 
        { 
            return it1->second;
        }
    else return *(new TermList ());
};

//
// Lattice
//

Lattice::Lattice():Terms(new TermStorage)
{
};

Lattice::~Lattice(){
delete Terms;
};

const Lattice::SiteMap& Lattice::getSiteMap() const
{
    return Sites;
}

const Lattice::TermStorage& Lattice::getTermStorage() const
{
    return *Terms;
}

void Lattice::printTerms(unsigned int n)
{
TermList Temp = Terms->getTerms(n);
for (TermList::const_iterator it1=Temp.begin(); it1!=Temp.end(); ++it1) {
    INFO(**it1 );
    };
}

void Lattice::printSites() const
{
    for (SiteMap::const_iterator it1=Sites.begin(); it1!=Sites.end(); ++it1) {
            INFO(*(it1->second));
        };
}

void Lattice::addSite(Lattice::Site* S)
{
    Sites[S->Label]= S ;
}

void Lattice::addSite(const std::string &Label, unsigned short orbitals, unsigned short spins)
{
    addSite(new Lattice::Site(Label, orbitals, spins));
}

void Lattice::addTerm(const Lattice::Term *T)
{
    unsigned int N=T->getOrder();
    for (unsigned int i=0; i<N; ++i) { // some checks to avoid bad terms.
        if (Sites.find(T->SiteLabels[i]) == Sites.end()) { ERROR("Sites, specified in this Term " << *T << " do not exist."); throw (exWrongLabel()); };
        if (T->Orbitals[i]>=Sites[T->SiteLabels[i]]->OrbitalSize) { ERROR("Wrong orbital indices in term " << *T << "."); throw (exWrongLabel()); };
        if (T->Spins[i]>=Sites[T->SiteLabels[i]]->SpinSize) { ERROR("Wrong spin indices in term " << *T << "."); throw (exWrongLabel()); };
    };
    if ( std::abs(T->Value) ) Terms->addTerm(T);
}

const Lattice::Site& Lattice::getSite(const std::string& Label) const
{
    std::map<std::string, Site*>::const_iterator it1=Sites.find(Label); 
    if (it1!=Sites.end()) throw (exWrongLabel()); 
    return *(it1->second);
}

const char* Lattice::exWrongLabel::what() const throw(){
    return "Wrong requested Label";
};

} // end of namespace Pomerol
