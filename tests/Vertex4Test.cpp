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

#include "Misc.h"
#include "Lattice.h"
#include "LatticePresets.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "OperatorPresets.h"
#include "IndexHamiltonian.h"
#include "Symmetrizer.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "FieldOperatorContainer.h"
#include "GFContainer.h"
#include "TwoParticleGFContainer.h"
#include "Vertex4.h"

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
    RealType omega3 = w(n3);

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

    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    U = 1.0;

    
    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&L, "A", U, -U/2.);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();
    IndexInfo.printIndices();

    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.diagonalize(world);

    srand(time(NULL));
    beta = 40.0 ;//+ 10.0*RealType(rand())/RAND_MAX;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, S, H);
    Operators.prepare();


    GreensFunction GF(S,H,Operators.getAnnihilationOperator(0), Operators.getCreationOperator(0), rho);
    GF.prepare();
    GF.compute();

    TwoParticleGF Chi_uuuu(S,H,Operators.getAnnihilationOperator(0), Operators.getAnnihilationOperator(0), Operators.getCreationOperator(0), 
                      Operators.getCreationOperator(0), rho);
    Chi_uuuu.prepare();
    Chi_uuuu.compute();
    
    ComplexType l,r;
    l = Chi_uuuu(0,0,0); r = gamma4ref_uuuu(0,0,0)*GF(0)*GF(0);
    if (!compare(l,r)) return EXIT_FAILURE;
    l = Chi_uuuu(2,5,2); r = gamma4ref_uuuu(2,5,2)*GF(2)*GF(5)*GF(5)*GF(2) - beta*GF(2)*GF(5);
    if (!compare(l,r)) return EXIT_FAILURE;
    l = Chi_uuuu(-10,-9,-10); r = gamma4ref_uuuu(-10,-9,-10)*GF(-10)*GF(-9)*GF(-10)*GF(-9) - beta*GF(-10)*GF(-9);
    if (!compare(l,r)) return EXIT_FAILURE;
    INFO("PASSED SUSC TEST");
    
    Vertex4 Gamma4_uuuu(Chi_uuuu,GF,GF,GF,GF);
    l = Gamma4_uuuu.value(2,5,2); r = gamma4ref_uuuu(2,5,2)*GF(2)*GF(5)*GF(5)*GF(2);
    if (!compare(l,r)) return EXIT_FAILURE;
    l = Gamma4_uuuu.value(-10,-9,-10); r = gamma4ref_uuuu(-10,-9,-10)*GF(-10)*GF(-9)*GF(-9)*GF(-10);
    if (!compare(l,r)) return EXIT_FAILURE;

    bool success = true; 
     for(int n1 = -10; n1<10; ++n1)
     for(int n2 = -10; n2<10; ++n2)
     for(int n3 = -10; n3<10; ++n3){
         l = Gamma4_uuuu.value(n1,n2,n3);
         r = gamma4ref_uuuu(n1,n2,n3)*GF(n1)*GF(n2)*GF(n3)*GF(n1+n2-n3);
         INFO(n1 << " " << n2 << " " << n3 << " " << n1+n2-n3 << " : " << l << " == " << r);
         success = success && compare(l,r);
     }
    if (!success) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
