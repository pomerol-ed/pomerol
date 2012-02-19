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
#include "LatticeAnalysis.h"
#include "IndexClassification.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "GFContainer.h"
#include "TwoParticleGFContainer.h"
#include "Vertex4Container.h"

#include<cstdlib>

using namespace Pomerol;

RealType U = 1.0;
RealType beta;

bool compare(ComplexType a, ComplexType b)
{
    return abs(a-b) < 1e-10;
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
            (1 + sqr(0.5*U/omega1))*
            (1 + sqr(0.5*U/omega2));

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
    Log.setDebugging(true);
    LatticeAnalysis Lattice;

    IndexClassification IndexInfo(Lattice);
    StatesClassification S(IndexInfo); 

    Hamiltonian H(IndexInfo,S);

    std::string LatticeFileName("gamma4.json");

    Lattice.readin(LatticeFileName);
    IndexInfo.prepare();
    S.compute();

    H.prepare();
    H.diagonalize();

    srand(time(NULL));
    beta = 40.0 ;//+ 10.0*RealType(rand())/RAND_MAX;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(S,H,IndexInfo);

    GFContainer G(IndexInfo,S,H,rho,Operators);

    std::set<IndexCombination2> GFindices;
    GFindices.insert(IndexCombination2(0,0));
    GFindices.insert(IndexCombination2(0,1));
    GFindices.insert(IndexCombination2(1,0));
    GFindices.insert(IndexCombination2(1,1));

    G.prepareAll(GFindices);
    G.computeAll(30);

    std::set<IndexCombination4> GF2indices;
    for(int i1=0; i1<=1; ++i1)
    for(int i2=0; i2<=1; ++i2)
    for(int i3=0; i3<=1; ++i3)
    for(int i4=0; i4<=1; ++i4)
        GF2indices.insert(IndexCombination4(i1,i2,i3,i4));

    TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
    Chi4.prepareAll(GF2indices);
    Chi4.computeAll(7);
    Vertex4Container Gamma4(Chi4,G);

    std::cout << Gamma4(0,0,0,0)(3,2,0) << std::endl;
//     std::cout << Chi4(0,0,0,0)(2,5,2) << std::endl;
//     std::cout << Chi4(0,0,0,0)(5,2,5) << std::endl;
//     std::cout << Chi4(0,0,0,0)(5,2,2) << std::endl;
//     std::cout << Chi4(0,0,0,0)(1,7,1) << std::endl;
//     std::cout << Chi4(0,0,0,0)(2,-2,4) << std::endl;
//     std::cout << Chi4(0,0,0,0)(29,-29,29) << std::endl;

//     std::cout << Chi4(0,1,0,1)(3,2,0) << std::endl;
//     std::cout << Chi4(0,1,0,1)(2,5,2) << std::endl;
//     std::cout << Chi4(0,1,0,1)(5,2,5) << std::endl;
//     std::cout << Chi4(0,1,0,1)(5,2,2) << std::endl;
//     std::cout << Chi4(0,1,0,1)(1,7,1) << std::endl;
//     std::cout << Chi4(0,1,0,1)(2,-2,4) << std::endl;
//     std::cout << Chi4(0,1,0,1)(29,-29,29) << std::endl;

//     Vertex4 Gamma4(IndexInfo,Chi4,G);
//     Gamma4.prepareUnAmputated();
//     Gamma4.computeUnAmputated();
//     Gamma4.prepareAmputated(GF2indices);
//     Gamma4.computeAmputated();
// 
//     for(int n1 = -10; n1<10; ++n1)
//     for(int n2 = -10; n2<10; ++n2)
//     for(int n3 = -10; n3<10; ++n3){
//         if( !compare(Gamma4.getAmputatedValue(IC(0,0,0,0),n1,n2,n3),gamma4ref_uuuu(n1,n2,n3)) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,1,1,1),n1,n2,n3),gamma4ref_uuuu(n1,n2,n3)) ||
//             !compare(Gamma4(IC(0,1,0,1),n1,n2,n3),gamma4ref_udud(n1,n2,n3)) ||
//             !compare(Gamma4(IC(1,0,1,0),n1,n2,n3),gamma4ref_udud(n1,n2,n3)) ||
//             !compare(Gamma4(IC(0,1,1,0),n1,n2,n3),-gamma4ref_udud(n2,n1,n3)) ||
//             !compare(Gamma4(IC(1,0,0,1),n1,n2,n3),-gamma4ref_udud(n2,n1,n3)) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,1,0,0),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(0,0,1,1),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,0,0,0),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(0,1,0,0),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(0,0,1,0),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(0,0,0,1),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(0,1,1,1),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,0,1,1),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,1,0,1),n1,n2,n3),0) ||
//             !compare(Gamma4.getAmputatedValue(IC(1,1,1,0),n1,n2,n3),0))
//         return EXIT_FAILURE;
//     }

    return EXIT_SUCCESS;
}
