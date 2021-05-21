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

using namespace Pomerol;

RealType U = 1.0;
RealType mu = 0.4;

// Reference Green's function
ComplexType Gref(int n, RealType beta)
{
    RealType omega = M_PI*(2*n+1)/beta;

    RealType w0 = 1.0;
    RealType w1 = exp(beta*mu);
    RealType w2 = exp(-beta*(-2*mu+U));
    RealType Z = w0 + 2*w1 +w2;
    w0 /= Z; w1 /= Z; w2 /= Z;

    return (w0+w1)/(I*omega+mu) + (w1+w2)/(I*omega+mu-U);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu);

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    IndexInfo.printIndices();

    INFO("Hamiltonian");
    INFO(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        INFO("Energy levels " << H.getEigenValues());
        INFO("The value of ground energy is " << H.getGroundEnergy());
    }

    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(rho.getWeight(i));

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex down_index = IndexInfo.getIndex("A",0,down);

    auto c_map = Operators.getCreationOperator(down_index).getBlockMapping();
    for (auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++)
    {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(down_index),
                      Operators.getCreationOperator(down_index),
                      rho);

    GF.prepare();
    GF.compute();

    auto compare = [](ComplexType a, ComplexType b) {
        return std::abs(a - b) < 1e-14;
    };

    bool result = true;
    for(int n = 0; n < 100; ++n) {
        INFO(GF(n) << " == " << Gref(n,beta));
        result = result && compare(GF(n), Gref(n,beta));
    }

    MPI_Finalize();

    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
