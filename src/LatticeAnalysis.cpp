#include "LatticeAnalysis.h"
#include <fstream>

std::ostream& operator<<(std::ostream& output, const SiteHoppingElement& out)
{
	output << out.From << "_{" << out.OrbitalFrom << "} -> " << out.To << "_{" << out.OrbitalTo << "} : " << out.Value;
	return output;
}
sLatticeSite::sLatticeSite (unsigned short type_, RealType filling_, unsigned short number_, RealType U):U(U)
{
	type = type_;
	number=number_;
	filling=filling_;
};

pLatticeSite::pLatticeSite(unsigned short type_, RealType filling_, unsigned short number_, RealType U, RealType J, string &basis):U(U),J(J),basis(basis)
{
	type = type_;
	number=number_;
	filling=filling_;
};

std::ostream& operator<<(std::ostream& output,const sLatticeSite& out)
{
	output << "Site N " << out.number << " is an s-orbital, filled by " << out.filling << " electrons, U = " << out.U << endl << "Hopping: " << endl;
	std::list<SiteHoppingElement*>::const_iterator it;
	for (it=out.HoppingList.begin();it!=out.HoppingList.end();++it){
		output << (**it) << endl;	
	}
return output;
}

std::ostream& operator<<(std::ostream& output,const pLatticeSite& out)
{
	output << "Site N " << out.number << " is a  p-orbital, filled by " << out.filling << " electrons, U = " << out.U << ", J = " << out.J << " in a " << out.basis << " basis"; 
	output << endl << "Hopping: " << endl;
	std::list<SiteHoppingElement*>::const_iterator it;
	for (it=out.HoppingList.begin();it!=out.HoppingList.end();++it){
		output << (**it) << endl;	
	}
return output;
}


LatticeAnalysis::LatticeAnalysis ()
{
	root = new Json::Value;
  	mapOrbitalValue["s"] = s;
  	mapOrbitalValue["p"] = p;
  	mapOrbitalValue["d"] = d;
  	mapOrbitalValue["f"] = f;
};

int LatticeAnalysis::readin()
{
  Json::Reader reader;
  std::ifstream in("Lattice.json");
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
		cout << ErrorException.what() << endl;
		exit(1);
  	}
  (*this).classifySites();
  return 0;
};

inline void LatticeAnalysis::enterHoppingListForCurrentSite(unsigned short CurrentSite, Json::Value &Hopping, std::list<SiteHoppingElement*> &HoppingList)
{
	for (unsigned short hop=0; hop<Hopping.size(); hop++){
		Json::Value CurrentHopping=Hopping[hop];
		unsigned short orbital_from = CurrentHopping["orbital_from"].asInt();
		unsigned short to = atoi(CurrentHopping["to"].asString().c_str());
		unsigned short orbital_to = CurrentHopping["orbital_to"].asInt();
		RealType t = CurrentHopping["value"].asDouble();
		SiteHoppingElement *H = new SiteHoppingElement(CurrentSite,orbital_from,to,orbital_to,t);
		HoppingList.push_back(H);
 		};
};


void LatticeAnalysis::classifySites()
{

  Json::Value sites = (*root)["sites"];
  for (unsigned short site = 0; site < sites.size(); site++ )
    {
	    std::stringstream current_site; current_site << site; 
	    switch (mapOrbitalValue[sites[current_site.str()]["type"].asString()])
	    {
		    case s: 
			    {
				    RealType U=sites[current_site.str()]["U"].asDouble();
					RealType filling = sites[current_site.str()]["filling"].asDouble();
					LatticeSite *S = new sLatticeSite(s,filling,site,U); 
  					Json::Value hopping = sites[current_site.str()]["hopping"];
					enterHoppingListForCurrentSite(site,hopping,S->HoppingList);
					SitesList.push_back(S);
				    break;
		            }

		    case p: 
			    {
				    RealType U=sites[current_site.str()]["U"].asDouble();
				    RealType J=sites[current_site.str()]["J"].asDouble();
					RealType filling = sites[current_site.str()]["filling"].asDouble();
				    string basis = sites[current_site.str()]["basis"].asString();
					LatticeSite *P = new pLatticeSite(p,filling,site,U,J,basis); 
  					Json::Value hopping = sites[current_site.str()]["hopping"];
					enterHoppingListForCurrentSite(site,hopping,P->HoppingList);
					SitesList.push_back(P);
				    break;
			    }

		    case d: break;
		    case f: break;
	    }

	};
};

const std::vector<LatticeSite*>& LatticeAnalysis::getSitesList()
{
	return SitesList;
}

std::stringstream& LatticeAnalysis::printSitesList()
{ 
	static std::stringstream current;
	std::vector<LatticeSite*>::const_iterator it;
	for (it=SitesList.begin();it!=SitesList.end();++it)
	{
			switch ((*it)->type)
			{
				case s: { current << (sLatticeSite&)(**it) << endl; break; };
				case p: { current << (pLatticeSite&)(**it) << endl; break; };
				default: { current << "ERROR"; break; };
			}
	}
	return current;
};
