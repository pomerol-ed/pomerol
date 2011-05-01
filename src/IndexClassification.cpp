// This file is part of pomerol ED code
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


#include "IndexClassification.h"

std::ostream& operator<<(std::ostream& output,const sIndexInfo& out)
{
output << "Index " << out.bitNumber << " of s-orbital, site N " << out.site << ", spin " << ((out.spin==1)?"up  ":"down") << ", U= " << out.U; 
return output;
}

std::ostream& operator<<(std::ostream& output,const pIndexInfo& out)
{
output << "Index " << out.bitNumber << " of p-orbital, site N " << out.site << ", spin " << ((out.spin==1)?"up  ":"down") << ", U= " << out.U << ", J=" << out.J; 
return output;
}

int IndexClassification::prepare()
{
  (*this).defineIndices();
  (*this).defineHopping();
  (*this).defineTerms();
  (*this).defineEquivalentIndexPermutations();

  return 0;
}


void IndexInfo::setIndexNumber(const unsigned short& in)
{
  bitNumber = in;
}

IndexClassification::IndexClassification(LatticeAnalysis &Lattice):Lattice(Lattice)
{
  N_bit=0;
}

void IndexClassification::defineIndices()
{
  INFO("IndexClassification: Defining Particle Indices according to the Lattice");
  std::vector<LatticeSite*> SitesList=Lattice.getSitesList();
  std::vector<LatticeSite*>::const_iterator SiteIterator;
  for (SiteIterator = SitesList.begin(); SiteIterator < SitesList.end(); SiteIterator++ )
    {
        switch ((*SiteIterator)->type)
        {
            case s: 
                {
                    sLatticeSite S = (sLatticeSite&)(**SiteIterator);
                    sIndexInfo *S1 = new sIndexInfo(S.number,s,0,S.U,S.LocalMu);
                    sIndexInfo *S2 = new sIndexInfo(S.number,s,1,S.U,S.LocalMu);
                    IndexInfoList.insert(IndexInfoList.begin()+N_bit/2,S1);
                    IndexInfoList.push_back(S2);

                    N_bit+=2; 
                    break;
                    }

            case p: 
                {
                    pLatticeSite P = (pLatticeSite&)(**SiteIterator);
                    pIndexInfo *P1 = new pIndexInfo(P.number,p,0,0,P.basis,P.U,P.J);
                    pIndexInfo *P2 = new pIndexInfo(P.number,p,0,1,P.basis,P.U,P.J);
                    pIndexInfo *P3 = new pIndexInfo(P.number,p,0,2,P.basis,P.U,P.J);
                    pIndexInfo *P4 = new pIndexInfo(P.number,p,1,0,P.basis,P.U,P.J);
                    pIndexInfo *P5 = new pIndexInfo(P.number,p,1,1,P.basis,P.U,P.J);
                    pIndexInfo *P6 = new pIndexInfo(P.number,p,1,2,P.basis,P.U,P.J);
                    IndexInfoList.insert(IndexInfoList.begin()+N_bit/2,P3);
                    IndexInfoList.insert(IndexInfoList.begin()+N_bit/2,P2);
                    IndexInfoList.insert(IndexInfoList.begin()+N_bit/2,P1);
                    IndexInfoList.push_back(P4);
                    IndexInfoList.push_back(P5);
                    IndexInfoList.push_back(P6);
                    N_bit+=6; 
                    break;
                }

            case d: N_bit+=10; break;
            case f: N_bit+=14; break;
        }

    }
  for (unsigned short bit = 0; bit < IndexInfoList.size(); bit++ ) IndexInfoList[bit]->setIndexNumber(bit);
}

void IndexClassification::defineTerms()
{
    INFO("IndexClassification: Defining all local terms");
   for (unsigned short bit = 0; bit < IndexInfoList.size()/2; bit++ ) 
      {
     switch (IndexInfoList[bit]->type)
        {
                case s:
             {
               sIndexInfo *current = (sIndexInfo*) IndexInfoList[bit];
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
             pIndexInfo *list[6];
             list[0]=(pIndexInfo*) IndexInfoList[bit];
             list[1]=(pIndexInfo*) IndexInfoList[bit+1];
             list[2]=(pIndexInfo*) IndexInfoList[bit+2];
             list[3]=(pIndexInfo*) IndexInfoList[bit+N_bit/2];
             list[4]=(pIndexInfo*) IndexInfoList[bit+N_bit/2+1];
             list[5]=(pIndexInfo*) IndexInfoList[bit+N_bit/2+2];

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

void IndexClassification::defineEquivalentIndexPermutations()
{
    INFO("IndexClassification: Defining equivalent ParticleIndex permutations");
    #warning Assuming explicitly spin symmetry - must be an input flag here 
    std::vector<ParticleIndex> CurrentAcceptedPermutation(N_bit); 
    // Spin symmetry
    for (int i=0;i<N_bit/2;i++) CurrentAcceptedPermutation[i]=i+N_bit/2;
    for (int i=N_bit/2;i<N_bit;i++) CurrentAcceptedPermutation[i]=i-N_bit/2;
    EquivalentPermutations.push_back(IndexPermutation(CurrentAcceptedPermutation));
}

void IndexClassification::printIndexInfoList()
{
  for (std::vector<IndexInfo*>::iterator it=IndexInfoList.begin(); it != IndexInfoList.end(); it++ ) (*it)->print_to_screen();
}

unsigned short IndexClassification::findIndex(const unsigned short &site, const unsigned short &orbital, const unsigned short &spin)
{
  for (unsigned short bit = 0; bit < IndexInfoList.size(); bit++ ) 
      if (IndexInfoList[bit]->site == site && IndexInfoList[bit]->spin == spin) 
        { 
          if (IndexInfoList[bit]->type == s) return bit;
          if (IndexInfoList[bit]->type == p)
            {
          pIndexInfo *current = (pIndexInfo*) IndexInfoList[bit];
          if ( ((current->basis == "Spherical") || (current->basis == "spherical" )) && ( current->orbital == orbital)) return bit;
        };
          /*
          if (IndexInfoList[bit]->type == "d")
            {
          if ( ( IndexInfoList[bit]->basis = "Spherical" || IndexInfoList[bit]->basis = "spherical" ) && ( IndexInfoList[bit]->orbital == 3)) return bit;
        };
          if (IndexInfoList[bit]->type == "f")
            {
          if ( ( IndexInfoList[bit]->basis = "Spherical" || IndexInfoList[bit]->basis = "spherical" ) && ( IndexInfoList[bit]->orbital == 5)) return bit;
        };
        */
        }
  return -1;
}

void IndexClassification::printHoppingMatrix()
{
    std::cout << HoppingMatrix << std::endl;
}

void IndexClassification::printTerms()
{
    std::cout << Terms << std::endl;
}

void IndexClassification::printEquivalentPermutations()
{
    INFO_NONEWLINE("(");
    for (ParticleIndex i=0; i<N_bit; ++i) INFO_NONEWLINE(i);
    INFO_NONEWLINE(")");
    INFO("");
    for (ParticleIndex i=0; i<N_bit+2; ++i) INFO_NONEWLINE("-");
    INFO("");
    for (std::list<IndexClassification::IndexPermutation>::const_iterator it=EquivalentPermutations.begin(); it!=EquivalentPermutations.end(); ++it) INFO(*it);
}
void IndexClassification::defineHopping()
/*
 * This procedure defines a Hopping Matrix from a Lattice input file
 */
{
    INFO("IndexClassification: Defining hopping terms");
    HoppingMatrix.resize(N_bit, N_bit);
    HoppingMatrix.setZero();
    std::vector<LatticeSite*> SitesList=Lattice.getSitesList();
    std::vector<LatticeSite*>::const_iterator SiteIterator;
    for (SiteIterator = SitesList.begin(); SiteIterator < SitesList.end(); SiteIterator++ ){
        std::list<SiteHoppingElement*> HoppingList = (*SiteIterator)->HoppingList; 
        std::list<SiteHoppingElement*>::const_iterator HoppingIterator;
        for (HoppingIterator=HoppingList.begin();HoppingIterator!=HoppingList.end();++HoppingIterator){
            for (unsigned short spin = 0; spin < 2; spin++){
                unsigned short bit_from = findIndex((*SiteIterator)->number, (*HoppingIterator)->OrbitalFrom, spin);
                unsigned short bit_to   = findIndex((*HoppingIterator)->To , (*HoppingIterator)->OrbitalTo  , spin);
                  HoppingMatrix(bit_from,bit_to)+=(*HoppingIterator)->Value;
            }
        }
    }
}

RealMatrixType& IndexClassification::getHoppingMatrix()
{
  return HoppingMatrix;
}

std::vector<IndexInfo*>& IndexClassification::getIndexInfoList()
{
  return IndexInfoList;
}
const int& IndexClassification::getIndexSize() const
{
    return N_bit;
}

bool IndexClassification::checkIndex(ParticleIndex in)
{
    return (in <= N_bit-1);
}

TermsList& IndexClassification::getTermsList()
{
    return Terms;
}

const std::list<IndexClassification::IndexPermutation>&  IndexClassification::getEquivalentIndexPermutations()
{
    return EquivalentPermutations;
}

void IndexClassification::definePorbitalSphericalTerms(pIndexInfo **list)
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
void IndexClassification::definePorbitalNativeTerms(pIndexInfo **list)
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
         output << **it << std::endl;
  }
return output;
}

IndexClassification::IndexPermutation::IndexPermutation(std::vector<ParticleIndex> &in)
{
    Permutation = in;
}

bool IndexClassification::IndexPermutation::isCorrect()
{
    bool correct = true;
    for (ParticleIndex i=0; i<Permutation.size() && correct; ++i)
                correct=(Permutation[i]<Permutation.size());
    for (ParticleIndex i=0; i<Permutation.size() && correct; ++i)
        for (ParticleIndex j=i+1; j<Permutation.size() && correct; ++j){
                correct=(Permutation[i]!=Permutation[j]);
            };
    return correct;
}

IndexClassification::IndexPermutation::~IndexPermutation()
{
    Permutation.clear();
}

std::ostream& operator<<(std::ostream& output, const IndexClassification::IndexPermutation& out)
{
    output << "(" << std::flush;
    for (ParticleIndex i=0; i<out.Permutation.size(); ++i) output << out.Permutation[i] << std::flush;
    output << ")" << std::flush;
    return output;
}
