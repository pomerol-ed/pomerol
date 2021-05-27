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
#include "LatticePresets.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operators.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"

using namespace Pomerol;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", 1.0, -0.5);
    HExpr += CoulombS("B", 2.0, -1.0);
    HExpr += Hopping("A", "B", -1.0);

    auto IndexInfo = MakeIndexClassification(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    INFO(H.getPart(BlockNumber(4)));
    H.compute(MPI_COMM_WORLD);
    INFO(H.getPart(BlockNumber(4)));
    RealType E = -2.8860009;
    RealType E_calc = H.getGroundEnergy();
    INFO("Lowest energy level is " << E_calc);
    if (std::abs(E-E_calc) > 1e-7) {
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

