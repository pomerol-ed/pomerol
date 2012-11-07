//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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

/** \file tests/IndexTerm.cpp
** \brief Test of the Operator::Term.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "Misc.h"
#include "Logger.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include <boost/shared_ptr.hpp>

using namespace Pomerol;

int main(int argc, char* argv[])
{
  /* Test of Operator::Term*/
  Log.setDebugging(true);
  bool Seq[4] = {1,0,1,0};
  ParticleIndex Ind[] = {0,0,1,1};
  std::vector<bool> seq_v(Seq, Seq+4);
  std::vector<ParticleIndex> ind_v(Ind, Ind+4);
  Operator::Term IT1(4, seq_v, ind_v, 1.0);
  INFO("Created Operator::Term" << IT1);
  INFO("Rearranging it to normal order");
  try {
    boost::shared_ptr<std::list<Operator::Term*> > out_terms=IT1.makeNormalOrder();
    INFO("Received " << IT1);
    INFO(out_terms->size() << " additional terms emerged : ");     
    for (std::list<Operator::Term*>::const_iterator it1=out_terms->begin(); it1!=out_terms->end(); ++it1) DEBUG(**it1);
    }
  catch (std::exception &e)
    {
        return EXIT_FAILURE;
    }
    

   FockState a1(4);
   a1[0]=1;
   a1[1]=0;

   seq_v.resize(1); ind_v.resize(1);
   seq_v[0]=1; ind_v[0]=1;
   Operator::Term IT2(1, seq_v, ind_v, 1.0);
   INFO("Acting with operator " << IT2 << " on a state " << a1 );
   FockState bra;
   MelemType result;
   boost::tie(bra, result) = IT2.actRight(a1);
   if (result != MelemType(-1)) return EXIT_FAILURE;
   else DEBUG ( "State: " << bra << " Matrix element: " << result);

   seq_v[0]=0; ind_v[0]=0;
   Operator::Term IT3 (1, seq_v, ind_v, 1.0);
   INFO("Acting with operator " << IT3 << " on a state " << a1 );
   boost::tie(bra, result) = IT3.actRight(a1);
   if (result != MelemType(1)) return EXIT_FAILURE;
   if (bra == ERROR_FOCK_STATE) DEBUG("Term vanishes")
   else DEBUG ( "State: " << bra << " Matrix element: " << result);

   seq_v.resize(2); ind_v.resize(2);
   seq_v[0]=1; seq_v[1]=0;
   ind_v[0]=1; ind_v[1]=1;
   Operator::Term IT4(2, seq_v, ind_v, 1.0);
   INFO("Checking term " << IT4);
   boost::shared_ptr<std::list<Operator::Term*> > out_terms=IT4.makeNormalOrder();
   INFO(out_terms->size() << " additional terms emerged : ");     
   for (std::list<Operator::Term*>::const_iterator it1=out_terms->begin(); it1!=out_terms->end(); ++it1) DEBUG(**it1);
   
   INFO("( " << IT4 << "==" << IT4 << " ) =" << (IT4==IT4));
   INFO("( " << IT4 << "==" << IT1 << " ) =" << (IT4==IT1));
   //DEBUG(IT4.commutes(IT4));

   seq_v[0]=1; seq_v[1]=0;
   ind_v[0]=0; ind_v[1]=1;
   Operator::Term IT5(2, seq_v, ind_v, 1.0);
   seq_v[0]=0; seq_v[1]=1;
   ind_v[0]=1; ind_v[1]=0;
   Operator::Term IT6(2, seq_v, ind_v, -1.0);
   Operator::Term IT7(2, seq_v, ind_v,  1.0);
   
   INFO("( " << IT5 << "==" << IT6 <<" ) = " << (IT5 == IT6));
   if (!(IT5 == IT6)) return EXIT_FAILURE;
   INFO("( " << IT5 << "==" << IT7 <<" ) = " << (IT5 == IT7));
   if ((IT5 == IT7)) return EXIT_FAILURE;
   
   seq_v.resize(4); ind_v.resize(4);
   seq_v[0]=1; seq_v[1]=0; seq_v[2]=1; seq_v[3]=0;
   ind_v[0]=0; ind_v[1]=1; ind_v[2]=2; ind_v[3]=3;
   Operator::Term IT8(4, seq_v, ind_v, 1.0);
   seq_v[0]=1; seq_v[1]=0; seq_v[2]=0; seq_v[3]=1;
   ind_v[0]=0; ind_v[1]=1; ind_v[2]=3; ind_v[3]=2;
   Operator::Term IT9(4, seq_v, ind_v, -1.0);
   INFO("( " << IT8 << "==" << IT9 <<" ) = " << (IT8 == IT9));
   if (!(IT8==IT9)) return EXIT_FAILURE;
   
   seq_v[0]=1; seq_v[1]=0; seq_v[2]=0; seq_v[3]=1;
   ind_v[0]=0; ind_v[1]=1; ind_v[2]=2; ind_v[3]=2;
   Operator::Term IT10(4, seq_v, ind_v, 1.0);
   seq_v[0]=0; seq_v[1]=1; seq_v[2]=0; seq_v[3]=1;
   ind_v[0]=1; ind_v[1]=0; ind_v[2]=2; ind_v[3]=2;
   Operator::Term IT11(4, seq_v, ind_v, -1.0);
   INFO("( " << IT10 << "==" << IT11 <<" ) = " << (IT10 == IT11));
   if (!(IT10==IT11)) return EXIT_FAILURE;
   
   INFO(IT4 << " commutes with " << IT4 << " = " << IT4.commutes(IT4));
   INFO(IT10 << " commutes with " << IT11 << " = " << IT10.commutes(IT11));

   // test reduce
   boost::shared_ptr<std::list<Operator::Term*> > list1(new std::list<Operator::Term*>);
   list1->push_back(new Operator::Term(IT10));
   list1->push_back(new Operator::Term(IT10));
   list1->push_back(new Operator::Term(IT11));
   list1->push_back(new Operator::Term(IT11));
   list1->push_back(new Operator::Term(IT10));
   list1->push_back(new Operator::Term(IT9));
   seq_v[0]=1; seq_v[1]=0; seq_v[2]=0; seq_v[3]=1;
   ind_v[0]=1; ind_v[1]=1; ind_v[2]=0; ind_v[3]=0;
   Operator::Term IT12 (4, seq_v, ind_v, 13.0);
   list1->push_back(new Operator::Term(IT12));
   Operator::Term IT12_2 (4, seq_v, ind_v, -5.0);
   Operator::Term IT12_3 (4, seq_v, ind_v, -8.0);
   list1->push_back(new Operator::Term(IT12_2));
   list1->push_back(new Operator::Term(IT12_3));
   INFO("Put 8 elements to list");

   Operator::Term::reduce(list1);
   INFO("Reduced to " << list1->size() << " elements.");
   if (list1->size()!=4) return EXIT_FAILURE;
   Operator::Term::prune(list1);
   INFO("Pruned to " << list1->size() << " elements.");
   if (list1->size()!=3) return EXIT_FAILURE;
   Operator A1(list1);
   A1.printAllTerms();

  /* end of test of Operator::Term */

  return EXIT_SUCCESS;
}

