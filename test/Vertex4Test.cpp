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

/** \file tests/gamma4.cpp
** \brief Test of a Green's function calculation (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
*/

#include "Misc.hpp"
#include "LatticePresets.hpp"
#include "Index.hpp"
#include "IndexClassification.hpp"
#include "Operators.hpp"
#include "StatesClassification.hpp"
#include "Hamiltonian.hpp"
#include "FieldOperatorContainer.hpp"
#include "GFContainer.hpp"
#include "TwoParticleGFContainer.hpp"
#include "Vertex4.hpp"

#include "./Utility.hpp"

#include<cstdlib>

using namespace Pomerol;

RealType U;
RealType beta;

bool compare(ComplexType a, ComplexType b, RealType tol = 1e-5)
{
    INFO("TEST: " << a << "?=" << b);
    return abs(a-b) < tol;
}

RealType delta(int n1, int n2)
{
    return (n1 == n2 ? 1.0 : 0.0);
}

RealType deltam(int n1, int n2)
{
    return (n1+n2==-1 ? 1.0 : 0.0);
}

inline RealType w(int n)
{
    return M_PI*(2*n+1)/beta;
}

// Reference gamma4(up,up,up,up)
ComplexType gamma4ref_uuuu(int n1, int n2, int n3)
{
    ComplexType Value = 0;

    RealType omega1 = w(n1);
    RealType omega2 = w(n2);

    Value = -beta*(delta(n1,n3)-delta(n2,n3))*sqr(0.5*U)*
            (1. + sqr(0.5*U/omega1))*
            (1. + sqr(0.5*U/omega2));

    return Value;
}

// Reference gamma4(up,down,up,down)
ComplexType gamma4ref_udud(int n1, int n2, int n3)
{
    ComplexType Value = 0;

    RealType omega1 = w(n1);
    RealType omega2 = w(n2);
    RealType omega3 = w(n3);
    RealType omega4 = omega1+omega2-omega3;

    RealType w = 1.0/(1.0+exp(beta*0.5*U));

    Value += U;
    Value += -0.125*U*U*U*(sqr(omega1)+sqr(omega2)+sqr(omega3)+sqr(omega4))/(omega1*omega2*omega3*omega4);
    Value += -0.1875*U*U*U*U*U/(omega1*omega2*omega3*omega4);
    Value += -beta*(2*deltam(n1,n2)+delta(n1,n3))*w*sqr(0.5*U)*(1.0+sqr(0.5*U/omega2))*(1.0+sqr(0.5*U/omega3));
    Value += beta*(2*delta(n2,n3)+delta(n1,n3))*(1-w)*sqr(0.5*U)*(1.0+sqr(0.5*U/omega1))*(1.0+sqr(0.5*U/omega2));

    return Value;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    U = 1.0;

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -U/2);

    auto IndexInfo = MakeIndexClassification(HExpr);
    std::cout << IndexInfo << std::endl;

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    srand(time(NULL));
    beta = 40.0 ;//+ 10.0*RealType(rand())/RAND_MAX;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();


    GreensFunction GF(S,H,Operators.getAnnihilationOperator(0), Operators.getCreationOperator(0), rho);
    GF.prepare();
    GF.compute();

    TwoParticleGF Chi_uuuu(S,H,Operators.getAnnihilationOperator(0), Operators.getAnnihilationOperator(0), Operators.getCreationOperator(0),
                      Operators.getCreationOperator(0), rho);
    Chi_uuuu.CoefficientTolerance = 1e-12;
    Chi_uuuu.prepare();
    Chi_uuuu.compute();

    ComplexType l,r;
    l = Chi_uuuu(0,0,0); r = gamma4ref_uuuu(0,0,0)*GF(0)*GF(0);
    if (!compare(l,r)) { MPI_Finalize(); return EXIT_FAILURE; }
    l = Chi_uuuu(2,5,2); r = gamma4ref_uuuu(2,5,2)*GF(2)*GF(5)*GF(5)*GF(2) - beta*GF(2)*GF(5);
    if (!compare(l,r)) { MPI_Finalize(); return EXIT_FAILURE; }
    l = Chi_uuuu(-10,-9,-10); r = gamma4ref_uuuu(-10,-9,-10)*GF(-10)*GF(-9)*GF(-10)*GF(-9) - beta*GF(-10)*GF(-9);
    if (!compare(l,r)) { MPI_Finalize(); return EXIT_FAILURE; }
    INFO("PASSED SUSC TEST");

    Vertex4 Gamma4_uuuu(Chi_uuuu,GF,GF,GF,GF);
    l = Gamma4_uuuu.value(2,5,2); r = gamma4ref_uuuu(2,5,2)*GF(2)*GF(5)*GF(5)*GF(2);
    if (!compare(l,r)) { MPI_Finalize(); return EXIT_FAILURE; }
    l = Gamma4_uuuu.value(-10,-9,-10); r = gamma4ref_uuuu(-10,-9,-10)*GF(-10)*GF(-9)*GF(-9)*GF(-10);
    if (!compare(l,r)) { MPI_Finalize(); return EXIT_FAILURE; }

    bool success = true;
     for(int n1 = -10; n1<10; ++n1)
     for(int n2 = -10; n2<10; ++n2)
     for(int n3 = -10; n3<10; ++n3){
         l = Gamma4_uuuu.value(n1,n2,n3);
         r = gamma4ref_uuuu(n1,n2,n3)*GF(n1)*GF(n2)*GF(n3)*GF(n1+n2-n3);
         INFO(n1 << " " << n2 << " " << n3 << " " << n1+n2-n3 << " : " << l << " == " << r);
         success = success && compare(l,r);
    }

    MPI_Finalize();
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
