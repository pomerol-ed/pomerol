#include "BitClassification.h"
#include <fstream>
#include <iomanip>

std::ostream& operator<<(std::ostream& output,const sBitInfo& out)
{
output << "Bit " << out.bitNumber << " of s-orbital, site N " << out.site << ", spin " << ((out.spin==1)?"up  ":"down") << ", U= " << out.U; 
return output;
}

std::ostream& operator<<(std::ostream& output,const pBitInfo& out)
{
output << "Bit " << out.bitNumber << " of p-orbital, site N " << out.site << ", spin " << ((out.spin==1)?"up  ":"down") << ", U= " << out.U << ", J=" << out.J; 
return output;
}

void BitInfo::setBitNumber(const unsigned short& in)
{
  bitNumber = in;
}

BitClassification::BitClassification()
{
  N_bit=0;
  mapOrbitalValue["s"] = s;
  mapOrbitalValue["p"] = p;
  mapOrbitalValue["d"] = d;
  mapOrbitalValue["f"] = f;
}

void BitClassification::defineBits()
{
  Json::Value sites = root["sites"];
  string type = "s";
  for (unsigned short site = 0; site < sites.size(); site++ )
    {
	    std::stringstream current_site; current_site << site; 
	    switch (mapOrbitalValue[sites[current_site.str()]["type"].asString()])
	    {
		    case s: 
			    {

				    string type = "s";
				    RealType U=sites[current_site.str()]["U"].asDouble();
				    sBitInfo *S1 = new sBitInfo(site,type,0,U);
				    sBitInfo *S2 = new sBitInfo(site,type,1,U);
				    BitInfoList.insert(BitInfoList.begin()+N_bit/2,S1);
				    BitInfoList.push_back(S2);
				    N_bit+=2; 
				    break;
		            }

		    case p: 
			    {
			    	    string type = "p";
				    RealType U=sites[current_site.str()]["U"].asDouble();
				    RealType J=sites[current_site.str()]["J"].asDouble();
				    string basis = sites[current_site.str()]["basis"].asString();
				    pBitInfo *P1 = new pBitInfo(site,type,0,0,basis,U,J);
				    pBitInfo *P2 = new pBitInfo(site,type,0,1,basis,U,J);
				    pBitInfo *P3 = new pBitInfo(site,type,0,2,basis,U,J);
				    pBitInfo *P4 = new pBitInfo(site,type,1,0,basis,U,J);
				    pBitInfo *P5 = new pBitInfo(site,type,1,1,basis,U,J);
				    pBitInfo *P6 = new pBitInfo(site,type,1,2,basis,U,J);
				    BitInfoList.insert(BitInfoList.begin()+N_bit/2,P3);
				    BitInfoList.insert(BitInfoList.begin()+N_bit/2,P2);
				    BitInfoList.insert(BitInfoList.begin()+N_bit/2,P1);
				    BitInfoList.push_back(P4);
				    BitInfoList.push_back(P5);
				    BitInfoList.push_back(P6);
				    N_bit+=6; 
				    break;
			    }

		    case d: N_bit+=10; break;
		    case f: N_bit+=14; break;
	    }

    }
  for (unsigned short bit = 0; bit < BitInfoList.size(); bit++ ) BitInfoList[bit]->setBitNumber(bit);
}

void BitClassification::printBitInfoList()
{
  for (vector<BitInfo*>::iterator it=BitInfoList.begin(); it != BitInfoList.end(); it++ ) (*it)->print_to_screen();
}
int BitClassification::readin()
{

  Json::Reader reader;
  std::ifstream in("Lattice.json");
  bool parsingSuccessful = reader.parse( in, root );
  if ( !parsingSuccessful )
  {
      // report to the user the failure and their locations in the document.
	std::cout  << "Failed to parse configuration\n";
	std::cout << reader.getFormatedErrorMessages();
        return 1;
  }

  (*this).defineBits();
  (*this).defineHopping();
  return 0;
}

unsigned short BitClassification::findBit(const unsigned short &site, const unsigned short &spin)
{
  for (unsigned short bit = 0; bit < BitInfoList.size(); bit++ ) 
	  if (BitInfoList[bit]->site == site && BitInfoList[bit]->spin == spin) 
	    { 
	      if (BitInfoList[bit]->type == "s") return bit;
	      if (BitInfoList[bit]->type == "p")
	        {
		  pBitInfo *current = (pBitInfo*) BitInfoList[bit];
		  if ( ((current->basis == "Spherical") || (current->basis == "spherical" )) && ( current->index == 1)) return bit;
		};
	      /*
	      if (BitInfoList[bit]->type == "d")
	        {
		  if ( ( BitInfoList[bit]->basis = "Spherical" || BitInfoList[bit]->basis = "spherical" ) && ( BitInfoList[bit]->index == 3)) return bit;
		};
	      if (BitInfoList[bit]->type == "f")
	        {
		  if ( ( BitInfoList[bit]->basis = "Spherical" || BitInfoList[bit]->basis = "spherical" ) && ( BitInfoList[bit]->index == 5)) return bit;
		};
		*/
	    }
  return -1;
}

vector<unsigned short>& BitClassification::findBits(const unsigned short &site) // This method find all bits, which correspond to a given site
{
  static std::vector<unsigned short> result;
  for (unsigned short bit = 0; bit < BitInfoList.size(); bit++ ) 
	  if (BitInfoList[bit]->site == site) result.push_back(bit);
  return result;
}

void BitClassification::printHoppingMatrix()
{
	cout << HoppingMatrix << endl;
}

//void BitClassification::createRecipe()
//{

//}
void BitClassification::defineHopping()
/*
 * This procedure defines a Hopping Matrix from a Lattice input file
 */
{

Json::Value sites = root["sites"];
HoppingMatrix.resize(N_bit, N_bit);
HoppingMatrix.setZero();
for (unsigned short from = 0; from < N_bit; from++ )
  {
	  unsigned short site_from = BitInfoList[from]->site;
	  std::stringstream current_site; current_site << site_from;
	  unsigned short spin_from = BitInfoList[from]->spin;
	  string type_from = BitInfoList[from]->type;

	  switch (mapOrbitalValue[type_from])
          {
	    case s:
	    {
	       for (unsigned short hoppingIndex=0; hoppingIndex < sites[current_site.str()]["hopping"].size(); hoppingIndex+=2)
	       {
		  unsigned short site_to = std::atoi(sites[current_site.str()]["hopping"][hoppingIndex].asString().c_str());
		  unsigned short to = (*this).findBit(site_to,spin_from); // This preserves spin conservation
		  HoppingMatrix(from,to)=sites[current_site.str()]["hopping"][hoppingIndex+1].asDouble();
	       };
	       break;
	    }
	    case p:
	    {
		string basis_from = sites[current_site.str()]["basis"].asString();
		pBitInfo *current = (pBitInfo*) BitInfoList[from];
	        if (( basis_from == "spherical" || basis_from == "Spherical") && ( current->index == 1))
	        	for (unsigned short hoppingIndex=0; hoppingIndex < sites[current_site.str()]["hopping"].size(); hoppingIndex+=2)
	       		{
		 		unsigned short site_to = std::atoi(sites[current_site.str()]["hopping"][hoppingIndex].asString().c_str());
		  		unsigned short to = (*this).findBit(site_to,spin_from); // This preserves spin conservation
		  		HoppingMatrix(from,to)=sites[current_site.str()]["hopping"][hoppingIndex+1].asDouble();
	       		};

		break;
	    }
	    default:
            { cout << "Undefined site " << site_from << " type." << endl; break; }
	  }
  }
}

RealMatrixType& BitClassification::getHoppingMatrix()
{
  return HoppingMatrix;
}

std::vector<BitInfo*>& BitClassification::getBitInfoList()
{
  return BitInfoList;
}
const int& BitClassification::getBitSize() const
{
	return N_bit;
}

