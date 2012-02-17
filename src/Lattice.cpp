#include "Lattice.h"
#include "LatticePresets.h"
#include <fstream>
#include <algorithm>

namespace Pomerol{

//
// Lattice::Site
//

Lattice::Site::Site()
{
};

Lattice::Site::Site(std::string label, unsigned short OrbitalSize, unsigned short SpinSize):label(label), OrbitalSize(OrbitalSize), SpinSize(SpinSize)
{
};


std::ostream& operator<<(std::ostream& output, const Lattice::Site& out)
{
    output << "Site \"" << out.label << "\", " << out.OrbitalSize << " orbital" << ((out.OrbitalSize>1)?"s":"") << ", " << out.SpinSize << " spin" << ((out.SpinSize>1)?"s":"") << ".";
	return output;
}

//
// Lattice::Term
//

Lattice::Term::Term (unsigned int N):N(N)
{
    Order.resize(N);
    SiteLabels.resize(N);
    Spins.resize(N);
    Orbitals.resize(N);
    for (unsigned int i=0; i<N; ++i) { Order[i]=false; SiteLabels[i]=""; Spins[i]=0; Orbitals[i]=0.0; }; 
    Value=0.0; 
};


Lattice::Term::Term(unsigned int N, bool * Order_, RealType Value_, std::string * SiteLabels_, unsigned short * Orbitals_, unsigned short *Spins_):
N(N)
{
  Order.assign( Order_, Order_+N );
  SiteLabels.assign( SiteLabels_, SiteLabels_+N );
  Spins.assign( Spins_, Spins_+N );
  Orbitals.assign( Orbitals_, Orbitals_+N );
  Value=Value_;
}

Lattice::Term::Term (const Lattice::Term &in):N(in.N), Order(in.Order), SiteLabels(in.SiteLabels), Spins(in.Spins), Orbitals(in.Orbitals), Value(in.Value)
{
};
unsigned int Lattice::Term::getOrder() const { return N; };

std::ostream& operator<< (std::ostream& output, const Lattice::Term& out)
{   
    output << out.Value << "*"; 
    for (unsigned int i=0; i<out.N; ++i) output << ((out.Order[i])?"c^{+}":"c") << "_{" << out.SiteLabels[i] << "," << out.Orbitals[i] << "," << out.Spins[i] << "}" ; 
    return output; 
};

//
// Lattice::TermStorage
//

Lattice::TermStorage::TermStorage()
{
};

int Lattice::TermStorage::addTerm(const Lattice::Term *T)
{
    unsigned int N = T->getOrder();
    Terms[N].push_back(new Term(*T));
    return 0;
};

const Lattice::TermList &Lattice::TermStorage::getTermList (unsigned int N)
{
    if (Terms.find(N)!=Terms.end()) 
        { 
            return Terms[N];
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
}

//
// JSONLattice
//

JSONLattice::JSONLattice(){
};


int JSONLattice::readin(const std::string &filename)
{
  Json::Value *root = new Json::Value;
  Json::Reader reader;
  std::ifstream in;
  in.open(filename.c_str());
  try
  {
    bool parsingSuccessful = reader.parse( in, *root );
    if ( !parsingSuccessful )
  	{
		std::cout  << "Failed to parse configuration\n";
		std::cout << reader.getFormatedErrorMessages();
        return 1;
  	}
  }
  catch (std::exception ErrorException)
  	{
		std::cout << ErrorException.what() << std::endl;
		exit(1);
  	}
  in.close();
  readSites((*root)["Sites"]);
  delete root;
  return 0;
};

void JSONLattice::readSites(Json::Value &JSONSites)
{
    for (Json::Value::iterator it=JSONSites.begin(); it!=JSONSites.end(); ++it){
        std::string label = it.key().asString();
        unsigned short OrbitalSize = (*it)["OrbitalSize"].asInt();
        DEBUG(OrbitalSize);
        DEBUG((*it)["U"].asDouble());
        }
    exit(0);
};

} // end of namespace Pomerol
