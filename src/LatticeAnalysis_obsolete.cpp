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


#include "LatticeAnalysis.h"
#include <fstream>

namespace Pomerol{

std::ostream& operator<<(std::ostream& output, const SiteHoppingElement& out)
{
	output << out.From << "_{" << out.OrbitalFrom << "} -> " << out.To << "_{" << out.OrbitalTo << "} : " << out.Value;
	return output;
}
sLatticeSite::sLatticeSite (unsigned short type_, RealType LocalMu_, unsigned short number_, RealType U):U(U)
{
	type = type_;
	number=number_;
	LocalMu=LocalMu_;
};

pLatticeSite::pLatticeSite(unsigned short type_, RealType LocalMu_, unsigned short number_, RealType U, RealType J, std::string &basis):U(U),J(J),basis(basis)
{
	type = type_;
	number=number_;
	LocalMu=LocalMu_;
};

std::ostream& operator<<(std::ostream& output,const sLatticeSite& out)
{
	output << "Site N " << out.number << " is an s-orbital, filled by " << out.LocalMu << " electrons, U = " << out.U << std::endl << "Hopping: " << std::endl;
	std::list<SiteHoppingElement*>::const_iterator it;
	for (it=out.HoppingList.begin();it!=out.HoppingList.end();++it){
		output << (**it) << std::endl;	
	}
return output;
}

std::ostream& operator<<(std::ostream& output,const pLatticeSite& out)
{
	output << "Site N " << out.number << " is a  p-orbital, filled by " << out.LocalMu << " electrons, U = " << out.U << ", J = " << out.J << " in a " << out.basis << " basis"; 
	output << std::endl << "Hopping: " << std::endl;
	std::list<SiteHoppingElement*>::const_iterator it;
	for (it=out.HoppingList.begin();it!=out.HoppingList.end();++it){
		output << (**it) << std::endl;	
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

int LatticeAnalysis::readin(std::string &LatticeFile)
{
  Json::Reader reader;
  std::ifstream in(LatticeFile.c_str());
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
					RealType LocalMu = sites[current_site.str()]["LocalMu"].asDouble();
					LatticeSite *S = new sLatticeSite(s,LocalMu,site,U); 
  					Json::Value hopping = sites[current_site.str()]["hopping"];
					enterHoppingListForCurrentSite(site,hopping,S->HoppingList);
					SitesList.push_back(S);
				    break;
		            }

		    case p: 
			    {
				    RealType U=sites[current_site.str()]["U"].asDouble();
				    RealType J=sites[current_site.str()]["J"].asDouble();
					RealType LocalMu = sites[current_site.str()]["LocalMu"].asDouble();
				    std::string basis = sites[current_site.str()]["basis"].asString();
					LatticeSite *P = new pLatticeSite(p,LocalMu,site,U,J,basis); 
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
			  case s: { current << (sLatticeSite&)(**it) << std::endl; break; };
			  case p: { current << (pLatticeSite&)(**it) << std::endl; break; };
				default: { current << "ERROR"; break; };
			}
	}
	return current;
};

} // end of namespace Pomerol
