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
   std::pair<FockState, RealType> result = IT2.act(a1);
   if (result.second != -1) return EXIT_FAILURE;
   else DEBUG ( "State: " << result.first << " Matrix element: " << result.second);

   seq_v[0]=0; ind_v[0]=0;
   Operator::Term IT3 (1, seq_v, ind_v, 1.0);
   INFO("Acting with operator " << IT3 << " on a state " << a1 );
   result = IT3.act(a1);
   if (result.second != 1) return EXIT_FAILURE;
   if (result.first == ERROR_FOCK_STATE) DEBUG("Term vanishes")
   else DEBUG ( "State: " << result.first << " Matrix element: " << result.second);
  /* end of test of Operator::Term */

  return EXIT_SUCCESS;
}

