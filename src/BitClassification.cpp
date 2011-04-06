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

int BitClassification::prepare()
{
  (*this).defineBits();
  (*this).defineHopping();
  (*this).defineTerms();
  return 0;
}


void BitInfo::setBitNumber(const unsigned short& in)
{
  bitNumber = in;
}

BitClassification::BitClassification(LatticeAnalysis &Lattice):Lattice(Lattice)
{
  N_bit=0;
}

void BitClassification::defineBits()
{
  std::vector<LatticeSite*> SitesList=Lattice.getSitesList();
  std::vector<LatticeSite*>::const_iterator SiteIterator;
  for (SiteIterator = SitesList.begin(); SiteIterator < SitesList.end(); SiteIterator++ )
    {
        switch ((*SiteIterator)->type)
        {
            case s: 
                {
                    sLatticeSite S = (sLatticeSite&)(**SiteIterator);
                    sBitInfo *S1 = new sBitInfo(S.number,s,0,S.U,S.LocalMu);
                    sBitInfo *S2 = new sBitInfo(S.number,s,1,S.U,S.LocalMu);
                    BitInfoList.insert(BitInfoList.begin()+N_bit/2,S1);
                    BitInfoList.push_back(S2);

                    N_bit+=2; 
                    break;
                    }

            case p: 
                {
                    pLatticeSite P = (pLatticeSite&)(**SiteIterator);
                    pBitInfo *P1 = new pBitInfo(P.number,p,0,0,P.basis,P.U,P.J);
                    pBitInfo *P2 = new pBitInfo(P.number,p,0,1,P.basis,P.U,P.J);
                    pBitInfo *P3 = new pBitInfo(P.number,p,0,2,P.basis,P.U,P.J);
                    pBitInfo *P4 = new pBitInfo(P.number,p,1,0,P.basis,P.U,P.J);
                    pBitInfo *P5 = new pBitInfo(P.number,p,1,1,P.basis,P.U,P.J);
                    pBitInfo *P6 = new pBitInfo(P.number,p,1,2,P.basis,P.U,P.J);
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
     switch (BitInfoList[bit]->type)
        {
                case s:
             {
               sBitInfo *current = (sBitInfo*) BitInfoList[bit];
             Term *T1 = new nnTerm(bit,bit+N_bit/2,current->U);
             Terms.addTerm(T1);    
             Term *N1 = new nTerm(bit,-current->LocalMu);
             Term *N2 = new nTerm(bit+N_bit/2,-current->LocalMu);
             Terms.addTerm(N1);
             Terms.addTerm(N2);
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

      //       Term *N1 = new nTerm(bit          ,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Term *N2 = new nTerm(bit+1        ,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Term *N3 = new nTerm(bit+2        ,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Term *N4 = new nTerm(bit+N_bit/2  ,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Term *N5 = new nTerm(bit+N_bit/2+1,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Term *N6 = new nTerm(bit+N_bit/2+2,(list[0]->U)*(-2.5) + 5. * (list[0]->J) );
      //       Terms.addTerm(N1); Terms.addTerm(N2); Terms.addTerm(N3); Terms.addTerm(N4); Terms.addTerm(N5); Terms.addTerm(N6);
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

unsigned short BitClassification::findBit(const unsigned short &site, const unsigned short &orbital, const unsigned short &spin)
{
  for (unsigned short bit = 0; bit < BitInfoList.size(); bit++ ) 
      if (BitInfoList[bit]->site == site && BitInfoList[bit]->spin == spin) 
        { 
          if (BitInfoList[bit]->type == s) return bit;
          if (BitInfoList[bit]->type == p)
            {
          pBitInfo *current = (pBitInfo*) BitInfoList[bit];
          if ( ((current->basis == "Spherical") || (current->basis == "spherical" )) && ( current->orbital == orbital)) return bit;
        };
          /*
          if (BitInfoList[bit]->type == "d")
            {
          if ( ( BitInfoList[bit]->basis = "Spherical" || BitInfoList[bit]->basis = "spherical" ) && ( BitInfoList[bit]->orbital == 3)) return bit;
        };
          if (BitInfoList[bit]->type == "f")
            {
          if ( ( BitInfoList[bit]->basis = "Spherical" || BitInfoList[bit]->basis = "spherical" ) && ( BitInfoList[bit]->orbital == 5)) return bit;
        };
        */
        }
  return -1;
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

    HoppingMatrix.resize(N_bit, N_bit);
    HoppingMatrix.setZero();
    std::vector<LatticeSite*> SitesList=Lattice.getSitesList();
    std::vector<LatticeSite*>::const_iterator SiteIterator;
    for (SiteIterator = SitesList.begin(); SiteIterator < SitesList.end(); SiteIterator++ ){
        std::list<SiteHoppingElement*> HoppingList = (*SiteIterator)->HoppingList; 
        std::list<SiteHoppingElement*>::const_iterator HoppingIterator;
        for (HoppingIterator=HoppingList.begin();HoppingIterator!=HoppingList.end();++HoppingIterator){
            for (unsigned short spin = 0; spin < 2; spin++){
                unsigned short bit_from = findBit((*SiteIterator)->number, (*HoppingIterator)->OrbitalFrom, spin);
                unsigned short bit_to   = findBit((*HoppingIterator)->To , (*HoppingIterator)->OrbitalTo  , spin);
                  HoppingMatrix(bit_from,bit_to)+=(*HoppingIterator)->Value;
            }
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

bool BitClassification::checkIndex(ParticleIndex in)
{
    return (in <= N_bit-1);
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
