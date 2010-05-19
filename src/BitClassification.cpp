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
  (*this).defineTerms();
  return 0;
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

void BitClassification::defineTerms()
{
   for (unsigned short bit = 0; bit < BitInfoList.size()/2; bit++ ) 
      {
	 switch (mapOrbitalValue[BitInfoList[bit]->type])
		{
	     	   case s:
		 	{
		  	 sBitInfo *current = (sBitInfo*) BitInfoList[bit];
			 Term *T1 = new nnTerm(bit,bit+N_bit/2,current->U);
			 Terms.addTerm(T1);	
			 break;
		 	};
		   case p:
		 	{
			 pBitInfo *list[6];
			 list[0]=(pBitInfo*) BitInfoList[bit];
			 list[1]=(pBitInfo*) BitInfoList[bit+1];
			 list[2]=(pBitInfo*) BitInfoList[bit+2];
			 list[3]=(pBitInfo*) BitInfoList[bit+N_bit/2];
			 list[4]=(pBitInfo*) BitInfoList[bit+N_bit/2+1];
			 list[5]=(pBitInfo*) BitInfoList[bit+N_bit/2+2];

			 // The array of all bits for p-orbital is created. First 3 bits are one spin direction, second 3 bits correspond to opposite direction.
			 if ((list[0]->basis == "spherical") || (list[0]->basis=="Spherical")) definePorbitalSphericalTerms(list); 
			if ((list[0]->basis == "native") || (list[0]->basis=="Native")) definePorbitalNativeTerms(list);

			 bit+=2;
			 break;
		        };
		   default : {break;};
	    	};
      };
}

void BitClassification::printBitInfoList()
{
  for (vector<BitInfo*>::iterator it=BitInfoList.begin(); it != BitInfoList.end(); it++ ) (*it)->print_to_screen();
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

void BitClassification::printTerms()
{
	cout << Terms << endl;
}
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



TermsList& BitClassification::getTermsList()
{
	return Terms;
}

void BitClassification::definePorbitalSphericalTerms(pBitInfo **list)
{
	RealType F0 = list[0]->U - 4./3.* list[1]->J; 
	RealType F2 = list[0]->J*25/3.;
	Eigen::Matrix3i  W1,W2,W3;
	
	W1 << 1,-2, 1,   -2, 4,-2,   1,-2, 1;
	W2 << 0, 3, 6,    3, 0, 3,   6, 3, 0;
	W3 << 0,-3, 0,   -3, 0,-3,   0,-3, 0;

	for (unsigned short m=0;m<3;m++)
	  for (unsigned short sigma=0; sigma<2; sigma++)
	    for (unsigned short m1=0;m1<3;m1++)
	      for (unsigned short sigma1=0; sigma1<2; sigma1++)
	{
		unsigned short bit_m_sigma = list[(sigma*3)+m]->bitNumber;
		unsigned short bit_m1_sigma1 = list[(sigma1*3)+m1]->bitNumber;
		
		if (m!=m1 && sigma == sigma1) { Term *T1 = new nnTerm(bit_m_sigma, bit_m1_sigma1, (F0-F2/5.)/2.); Terms.addTerm(T1); }
		if (sigma != sigma1) { Term *T2 = new nnTerm(bit_m_sigma, bit_m1_sigma1, (F0)/2.); Terms.addTerm(T2); }
		if (sigma != sigma1) { Term *T3 = new nnTerm(bit_m_sigma, bit_m1_sigma1, (F2)/2./25.*W1(m,m1)); if (W1(m,m1)!=0) Terms.addTerm(T3); }

		if (sigma != sigma1)
		{
		unsigned short bit_m_sigma1 = list[(sigma1*3)+m]->bitNumber;
		unsigned short bit_m1_sigma = list[(sigma*3)+m1]->bitNumber;
		Term *T4 = new spinflipTerm(bit_m_sigma,bit_m1_sigma1, bit_m_sigma1, bit_m1_sigma, F2/2./25.*W2(m,m1));
		if (W2(m,m1)!=0) Terms.addTerm(T4);

		unsigned short bit_minusm_sigma1=bit_m_sigma1 = list[(sigma1*3)+2-m]->bitNumber;
		unsigned short bit_minusm1_sigma=bit_m_sigma1 = list[(sigma*3)+2-m1]->bitNumber;
		Term *T5 = new spinflipTerm(bit_m_sigma,bit_minusm_sigma1, bit_m1_sigma1, bit_minusm1_sigma, F2/2./25.*W3(m,m1));
		if (W3(m,m1)!=0) Terms.addTerm(T5);
		}
	}
};
void BitClassification::definePorbitalNativeTerms(pBitInfo **list)
{
	RealType U = list[0]->U;
	RealType J = list[0]->J;
	for (unsigned short p=0;p<3;p++)
	  for (unsigned short sigma=0; sigma<2; sigma++)
	    for (unsigned short p1=0;p1<3;p1++)
	{
		unsigned short sigma1 = 1-sigma;

		unsigned short bit_p_sigma = list[(sigma*3)+p]->bitNumber;
		unsigned short bit_p_sigma1 = list[(sigma1*3)+p]->bitNumber;
		Term *T = new nnTerm(bit_p_sigma,bit_p_sigma1,U/2.);
		Terms.addTerm(T);
		
		if (p!=p1)
		{
			unsigned short bit_p1_sigma1 = list[(sigma1*3)+p1]->bitNumber;
			Term *T2 = new nnTerm(bit_p_sigma,bit_p1_sigma1,(U-2.*J)/2.);
			Terms.addTerm(T2);

			unsigned short bit_p1_sigma = list[(sigma*3)+p1]->bitNumber;
			Term *T3 = new nnTerm(bit_p_sigma,bit_p1_sigma,(U-3.*J)/2.);
			Terms.addTerm(T3);

			Term *T4 = new spinflipTerm(bit_p_sigma,bit_p1_sigma1, bit_p1_sigma, bit_p_sigma1,(-1.*J)/2.);
			Terms.addTerm(T4);

			Term *T5 = new spinflipTerm(bit_p1_sigma,bit_p1_sigma1, bit_p_sigma, bit_p_sigma1,(-1.*J)/2.);
			Terms.addTerm(T5);
		}
	}
};


void TermsList::addTerm(Term* in)
{
	unsigned short order = in->N;
	TermsMap[in->N].push_back(in);
	maxOrder=(maxOrder>=order)?maxOrder:order;
}

std::list<Term*>& TermsList::getTerms(unsigned short order)
{
	return TermsMap[order];
}

std::ostream& operator<<(std::ostream& output, TermsList& out)
{
  for (unsigned short order=2;order<=out.maxOrder;order+=2)
  {
     std::list<Term*>::iterator it;
     for ( it=(out.getTerms(order)).begin() ; it != out.getTerms(order).end(); it++ ) 
	     output << **it << endl;
  }
return output;
}
