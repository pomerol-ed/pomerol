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
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "OperatorPresets.h"


using namespace Pomerol;
using namespace Pomerol::OperatorPresets;

int main(int argc, char* argv[])
{
  /* Test of Operator::Term*/

  Operator IT1 = Cdag(0)*C(0)*Cdag(1)*C(1);
  INFO("Created Operator " << IT1);

   FockState a1(4);
   a1[0]=1;
   a1[1]=0;
   FockState res_state;
   MelemType result;
   std::map<FockState, MelemType> out;

   Operator IT2 = Cdag(1);
   out=IT2.actRight(a1);
   res_state = out.begin()->first;
   result = out.begin()->second;
   INFO ( IT2 << "|" << a1 << "> =" << result << "|" << res_state << ">");
   if (result != MelemType(-1)) return EXIT_FAILURE;

   Operator IT3 = C(0);
   out=IT3.actRight(a1);
   res_state = out.begin()->first;
   result = out.begin()->second;
   INFO ( IT3 << "|" << a1 << "> =" << result << "|" << res_state << ">");
   if (result != MelemType(1)) return EXIT_FAILURE;
   if (res_state == ERROR_FOCK_STATE) DEBUG("Term vanishes");

   INFO(IT3 << "*" << IT2 << " = " << IT3*IT2);
   INFO(IT2 << "*" << IT3 << " = " << IT2*IT3);
   INFO("(" << IT2 << "*" << IT3 << "==" << IT3 << "*" << IT2 << " ) = " << (IT2*IT3 == IT3*IT2));
   INFO(IT2 << " commutes with " << IT3 << " = " << IT2.commutes(IT3));
   if (IT2.commutes(IT3)) return EXIT_FAILURE;

   Operator IT4 = Cdag(1)*C(1);
   INFO("( " << IT4 << "==" << IT4 << " ) = " << (IT4==IT4));
   INFO("( " << IT4 << "==" << IT1 << " ) = " << (IT4==IT1));
   if (IT1 == IT4) return EXIT_FAILURE;
   INFO(IT4 << " commutes with " << IT4 << " = " << IT4.commutes(IT4));
   if (!(IT4.commutes(IT4))) return EXIT_FAILURE;

   Operator IT5 = Cdag(0)*C(1);

   Operator IT6 = - C(1)*Cdag(0);
   Operator IT7 = C(1)*Cdag(0);

   INFO("( " << IT5 << "==" << IT6 <<" ) = " << (IT5 == IT6));
   if (!(IT5 == IT6)) return EXIT_FAILURE;
   INFO("( " << IT5 << "==" << IT7 <<" ) = " << (IT5 == IT7));
   if ((IT5 == IT7)) return EXIT_FAILURE;

   Operator IT7_2 = Cdag(2)*C(2);
   INFO(IT4 << " commutes with " << IT7_2 << " = " << IT4.commutes(IT7_2));

   Operator IT101 = Cdag(1)*C(1); //((std::make_tuple(1.0, ops)));

   Operator IT102 = Cdag(2)*C(2)*Cdag(0)*C(0); //((std::make_tuple(1.0, ops)));
   if (!(IT102.commutes(IT101))) return EXIT_FAILURE;

   Operator IT8 = Cdag(0)*C(1)*Cdag(2)*C(3);

   Operator IT9 = - Cdag(0) * C(1) * C(3) * Cdag(2); // (std::make_tuple(-1.0, ops));
   INFO("( " << IT8 << "==" << IT9 <<" ) = " << (IT8 == IT9));
   if (!(IT8==IT9)) return EXIT_FAILURE;
   //ops[0]=std::make_tuple(1,0);
   //ops[1]=std::make_tuple(0,1);
   //ops[2]=std::make_tuple(0,2);
   //ops[3]=std::make_tuple(1,2);
   Operator IT10 = Cdag(0)*C(1)*C(2)*Cdag(2);

   Operator IT11 = -C(1)*Cdag(0)*C(2)*Cdag(2);
   INFO("( " << IT10 << "==" << IT11 <<" ) = " << (IT10 == IT11));
   if (!(IT10==IT11)) return EXIT_FAILURE;

   INFO(IT10 << " commutes with " << IT11 << " = " << IT10.commutes(IT11));

   // test reduce
   Operator IT12Base = Cdag(1)*C(1)*C(0)*Cdag(0);
   Operator IT12 = 13*IT12Base; //(std::make_tuple(13.0, ops));
   Operator IT12_2 = -5 * IT12Base; // (std::make_tuple(-5.0, ops));
   Operator IT12_3 = -8*IT12Base; //(std::make_tuple(-8.0, ops));

   Operator ITsum;
   ITsum=IT12;
   ITsum+=IT12_2;
   ITsum+=IT12_3;

   INFO(IT12 << "+" << IT12_2 << "+" << IT12_3 << "=" << ITsum);
   if (!ITsum.isEmpty()) return EXIT_FAILURE;

   Operator Op1 = C(0)*Cdag(0)*C(1)*Cdag(1);
   Operator Op2 = Cdag(0)*C(0)*Cdag(1)*C(1);
   for (size_t i=0; i<=4; i++) {
       INFO("<" << FockState(4,i) << "|" << Op1 << "|" << FockState(4,i) << "> = " << Op1.getMatrixElement(FockState(4,i),FockState(4,i)));
       INFO("<" << FockState(4,i) << "|" << Op2 << "|" << FockState(4,i) << "> = " << Op2.getMatrixElement(FockState(4,i),FockState(4,i)));
        };

  /* end of test of Operator::Term */

  return EXIT_SUCCESS;
}

