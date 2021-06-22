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

/** \file tests/EnsembleAverage1siteTest.cpp
** \brief Test of a calculation of ensemble averages (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
*/

#include "Misc.hpp"
#include "LatticePresets.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "Operators.hpp"
#include "StatesClassification.hpp"
#include "Hamiltonian.hpp"
#include "FieldOperatorContainer.hpp"
#include "EnsembleAverage.hpp"

#include "./Utility.hpp"

#include <cmath>
#include <cstdlib>

using namespace Pomerol;

RealType U = 1.0;
RealType mu = 0.4;
RealType h_field = 0.01;
//RealType h_field = 0;

bool compare(ComplexType a, ComplexType b)
{
    return std::abs(a-b) < 1e-14;
}

// Exact result
struct ExactResult{
    RealType n_up, n_down;

    ExactResult(RealType beta){
        RealType w0 = 1.0;
        RealType wu = exp(beta*(mu+h_field));
        RealType wd = exp(beta*(mu-h_field));
        RealType w2 = exp(-beta*(-2*mu+U));
        RealType Z = w0 + wu + wd + w2;
        w0 /= Z; w2 /= Z; wu /= Z; wd /= Z;

        n_up = wu + w2;
        n_down = wd + w2;
    }
};

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    // h_field (n_down - n_up)
    // addHopping double-counts diagonal term, and divide h_field by 2.
    auto HExpr = CoulombS("A", U, -mu) + Magnetization("A", -h_field);

    INFO("Hamiltonian");
    INFO(HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    std::cout << IndexInfo << std::endl;

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex dn_index = IndexInfo.getIndex("A",0,down);
    ParticleIndex up_index = IndexInfo.getIndex("A",0,up);

    // quadratic operators, c^+ c
    QuadraticOperator s_plus(IndexInfo, HS, S, H, up_index, dn_index);
    QuadraticOperator s_minus(IndexInfo, HS, S, H, dn_index, up_index);
    QuadraticOperator n_up(IndexInfo, HS, S, H, up_index, up_index);
    QuadraticOperator n_dn(IndexInfo, HS, S, H, dn_index, dn_index);

    std::vector<QuadraticOperator*> quad_ops{&s_plus, &s_minus, &n_up, &n_dn};
    for(auto op : quad_ops){
        op->prepare(HS);
        op->compute();
    }

    // for print
    std::vector<std::string> names = {"< S_+ >", "< S_- >", "< n_up >", "< n_down >"};

    // reference data
    ExactResult Exact(beta);
    std::vector<RealType> Refs = {0, 0, Exact.n_up, Exact.n_down};

    // compute susceptibilities, and compare them with reference data
    bool result = true;
    for(int i=0; i<quad_ops.size(); i++){
        print_section(names[i]);
        EnsembleAverage EA(S,H, *quad_ops[i], rho);
        EA.compute();

        // check if results are correct
        INFO(EA.getResult() << " == " << Refs[i]);
        result = (result && compare(EA.getResult(), Refs[i]));
    }

    MPI_Finalize();
    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
