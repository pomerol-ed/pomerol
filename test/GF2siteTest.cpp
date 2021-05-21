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

/** \file tests/green.cpp
** \brief Test of a Green's function calculation (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
*/

#include "Misc.h"
#include "LatticePresets.h"
#include "Operators.h"
#include "IndexClassification.h"
#include "HilbertSpace.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "DensityMatrix.h"
#include "FieldOperatorContainer.h"
#include "GFContainer.h"

#include "./Utility.h"

#include <cstdlib>
#include <tuple>
#include <vector>

using namespace Pomerol;

RealType U = 1.0;
RealType mu = 0.5;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu);
    HExpr += CoulombS("B", U, -mu);
    HExpr += Hopping("A","B", -1.0);

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    IndexInfo.printIndices();

    auto A_dn = IndexInfo.getInfo(0);
    auto A_up = IndexInfo.getInfo(1);
    auto B_dn = IndexInfo.getInfo(2);
    auto B_up = IndexInfo.getInfo(3);

    auto N = Operators::N(std::vector<decltype(A_up)>{A_up, A_dn, B_up, B_dn});
    INFO("N = " << N);

    auto Sz = Operators::Sz(std::vector<decltype(A_up)>{A_up, B_up},
                            std::vector<decltype(A_dn)>{A_dn, B_dn});
    INFO("Sz = " << Sz);

    if(Sz * N != N * Sz) return EXIT_FAILURE;
    if(HExpr * N != N * HExpr) return EXIT_FAILURE;
    INFO("H commutes with N");
    if(HExpr * Sz != Sz * HExpr) return EXIT_FAILURE;
    INFO("H commutes with Sz");

    INFO("Hamiltonian");
    INFO(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.getPart(BlockNumber(4)).print_to_screen();
    H.getPart(BlockNumber(5)).print_to_screen();
    H.compute(MPI_COMM_WORLD);
    H.getPart(BlockNumber(4)).print_to_screen();
    H.getPart(BlockNumber(5)).print_to_screen();
    INFO("The value of ground energy is " << H.getGroundEnergy());

    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(rho.getWeight(i));

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    auto c_map = Operators.getCreationOperator(0).getBlockMapping();
    for(auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++)
    {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(0),
                      Operators.getCreationOperator(0),
                      rho);

    GF.prepare();
    GF.compute();

    ComplexVectorType G_ref(10);
    G_ref << -2.53021005e-01*I,
             -4.62090702e-01*I,
             -4.32482782e-01*I,
             -3.65598615e-01*I,
             -3.07785174e-01*I,
             -2.62894141e-01*I,
             -2.28274316e-01*I,
             -2.01170772e-01*I,
             -1.79539602e-01*I,
             -1.61950993e-01*I;

    bool result = true;
    auto compare = [](ComplexType a, ComplexType b) {
        // The tolerance has to be fairly large as Pomerol discards
        // some contributions to the GF.
        return std::abs(a - b) < 1e-8;
    };

    for(int n = 0; n<10; ++n) {
        INFO(GF(n) << " == " << G_ref(n));
        result = result && compare(GF(n), G_ref(n));
    }

    MPI_Finalize();

    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
