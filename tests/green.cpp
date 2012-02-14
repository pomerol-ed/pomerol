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

/** \file tests/green.h
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

#include<cstdlib>

using namespace Pomerol;

RealType U = 1.0;
RealType mu = 0.4;

bool compare(ComplexType a, ComplexType b)
{
    return abs(a-b) < 1e-14;
}

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
    LatticeAnalysis Lattice;

    IndexClassification IndexInfo(Lattice);
    StatesClassification S(IndexInfo); 

    Hamiltonian H(IndexInfo,S);

    std::string LatticeFileName("green.json");

    Lattice.readin(LatticeFileName);
    IndexInfo.prepare();
    S.compute();

    H.prepare();
    H.compute();

    srand (time(NULL));
    RealType beta = 10.0 + 10.0*RealType(rand())/RAND_MAX;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(S,H,IndexInfo);
    
    GFContainer G(S,H,rho,IndexInfo,Operators);

    std::vector<GFContainer::IndexCombination*> indices;
    indices.push_back(new GFContainer::IndexCombination(0,0));
    indices.push_back(new GFContainer::IndexCombination(0,1));
    indices.push_back(new GFContainer::IndexCombination(1,0));
    indices.push_back(new GFContainer::IndexCombination(1,1));
    
    G.readInitialIndices(indices);
    G.prepare();
    G.compute();

    for(int n = -100; n<100; ++n)
        if( !compare(G(0,0,n),Gref(n,beta)) ||
            !compare(G(0,1,n),0.0) ||
            !compare(G(1,0,n),0.0) ||
            !compare(G(1,1,n),Gref(n,beta))
            )
            return EXIT_FAILURE;

    return EXIT_SUCCESS;
}