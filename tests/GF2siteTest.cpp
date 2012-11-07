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
#include "Logger.h"
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
    return abs(a-b) < 1e-5;
}

// Reference Green's function

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}



int main(int argc, char* argv[])
{
    Log.setDebugging(true);
    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&L, "A", 1.0, -0.5);
    L.addSite(new Lattice::Site("B",1,2));
    LatticePresets::addCoulombS(&L, "B", 1.0, -0.5);

    LatticePresets::addHopping(&L, "A","B", -1.0);
    INFO("Sites");
    L.printSites();
    INFO("Terms with 2 operators");
    L.printTerms(2);
    INFO("Terms with 4 operators");
    L.printTerms(4);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();
    print_section("Indices");
    IndexInfo.printIndices();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();
    INFO("Terms with 2 operators");
    Storage.printTerms(2);
    INFO("Terms with 4 operators");
    Storage.printTerms(4);
    
    print_section("Making hamiltonian normal ordered");
    Storage.printAllTerms();
    INFO("==");
    Storage.makeNormalOrder();
    INFO("==");
    Storage.printAllTerms();


    ParticleIndex IndexSize = IndexInfo.getIndexSize();
    OperatorPresets::Sz Sz(IndexSize);
    INFO("Sz terms");
    Sz.printAllTerms();
    OperatorPresets::N N(IndexSize);
    INFO("N terms");
    N.printAllTerms();
    INFO(Sz.commutes(N));
    exit(0);

    DEBUG(N);
    if (!(Storage.commutes(N))) return EXIT_FAILURE;
    INFO("H commutes with N");
    if (!(Storage.commutes(Sz))) return EXIT_FAILURE;
    INFO("H commutes with Sz");
    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.getPart(BlockNumber(4)).print_to_screen();
    H.getPart(BlockNumber(5)).print_to_screen();
    H.diagonalize();
    H.getPart(BlockNumber(4)).print_to_screen();
    H.getPart(BlockNumber(5)).print_to_screen();
    INFO("The value of ground energy is " << H.getGroundEnergy());

    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(H.getEigenValue(i)); 
 
    //srand (time(NULL));
    RealType beta = 10.0; // + 10.0*RealType(rand())/RAND_MAX;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(rho.getWeight(i)); 

    FieldOperatorContainer Operators(IndexInfo, S, H);
    Operators.prepare();

    auto c_map=Operators.getCreationOperator(0).getNonTrivialIndices();
    for (auto c_map_it=c_map.begin(); c_map_it!=c_map.end(); c_map_it++)
        {
            INFO(c_map_it->second << "->" << c_map_it->first);
        }

    GreensFunction GF(S,H,Operators.getAnnihilationOperator(0), Operators.getCreationOperator(0), rho);

    GF.prepare();
    GF.compute(10);

    RealVectorType GF_im(10);
    GF_im<<-2.53021005e-01, -4.62090702e-01, -4.32482782e-01, -3.65598615e-01, -3.07785174e-01, -2.62894141e-01, 
           -2.28274316e-01, -2.01170772e-01, -1.79539602e-01, -1.61950993e-01;
 
    for(int n = 0; n<10; ++n) {
        DEBUG(GF(n) << " " << GF_im(n)*I);
        if( !compare(GF(n),GF_im(n)*I))
            return EXIT_FAILURE;
        }

//    for(int n = -100; n<100; ++n)
//        if( !compare(GF(n),Gref(n,beta)))
//            return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
