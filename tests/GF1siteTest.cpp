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

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    
    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));

    LatticePresets::addCoulombS(&L, "A", U, -mu);
    print_section("Sites");
    L.printSites();
    print_section("Terms");
    L.printTerms(2);
    INFO("Terms with 4 operators");
    L.printTerms(4);


    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();
    print_section("Indices");
    IndexInfo.printIndices();

    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.compute(world);
 
    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, S, H);
    Operators.prepare();
    Operators.compute();

    ParticleIndex down_index = IndexInfo.getIndex("A",0,down);
    
    FieldOperator::BlocksBimap c_map = Operators.getCreationOperator(down_index).getBlockMapping();
    for (FieldOperator::BlocksBimap::right_const_iterator c_map_it=c_map.right.begin(); c_map_it!=c_map.right.end(); c_map_it++)
        {
            INFO(c_map_it->first << "->" << c_map_it->second);
        }

    GreensFunction GF(S,H,Operators.getAnnihilationOperator(down_index), Operators.getCreationOperator(down_index), rho);

    GF.prepare();
    GF.compute();

    bool result = true;
    for(int n = 0; n<100; ++n) {
        INFO(GF(n) << " == " << Gref(n,beta));
        result = (result && compare(GF(n),Gref(n,beta)));
        }
    if (!result) return EXIT_FAILURE;
    return EXIT_SUCCESS;

    return EXIT_SUCCESS;
}
