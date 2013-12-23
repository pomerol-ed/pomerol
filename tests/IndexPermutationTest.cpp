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
** \brief Test of the Symmetrizer::IndexPermutation.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "IndexHamiltonian.h"
#include "Symmetrizer.h"
#include <boost/shared_ptr.hpp>

using namespace Pomerol;

int main(int argc, char* argv[])
{
  /* Make a check of the Symmetrizer::IndexPermutation. */
  ParticleIndex tmp[]={1,2,3,0,4,5};
  std::vector<ParticleIndex> tmp2(tmp, tmp+6);
  DynamicIndexCombination A(tmp2);
  /* end of test of Symmetrizer::IndexPermutation. */
  try {
        Symmetrizer::IndexPermutation B(A);
        INFO("Constructed a permutation of indices " << A );
        INFO("The cycle length of it is " << B.getCycleLength());
    }
  catch (std::exception &e)
    {
        return EXIT_FAILURE;
    }
  /* end of test of IndexHamiltonian::Term */

  return EXIT_SUCCESS;
}

