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

#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include <boost/filesystem/path.hpp>

#include<cstdlib>
#include <fstream>

#pragma clang diagnostic ignored "-Wc++11-extensions"

using namespace Pomerol;

bool compare(ComplexType a, ComplexType b)
{
    return abs(a-b) < 1e-5;
}

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-F1(y))<tolerance);
}



int main(int argc, char* argv[])
{
    Log.setDebugging(true);
    Lattice L;
    print_section("Kondo chain diagonalization");
    
    size_t NSites = 1;
    RealType t = 0.0;
    RealType U = 100;
    RealType mu = U/2;
    RealType J = 2.3;

    INFO("Diagonalization of " << NSites << "+1 sites");
    L.addSite(new Lattice::Site("K",1,2));
    std::map<size_t, std::string> chain_site_names;
    LatticePresets::addCoulombS(&L, "K", U, -U/2.);
    for (size_t i=0; i<NSites; ++i) {
        std::stringstream s;
        s << "T" << i;
        chain_site_names[i]=s.str();
        L.addSite(new Lattice::Site(chain_site_names[i],1,2));
        LatticePresets::addLevel(&L, chain_site_names[i], -mu);
        }

    INFO("Sites");
    L.printSites();

    for (size_t i=0; i<NSites-1; ++i) {
        LatticePresets::addHopping(&L, chain_site_names[i],chain_site_names[i+1], -t);
    }
    if (NSites>2) LatticePresets::addHopping(&L, chain_site_names[0],chain_site_names[NSites-1], -t);
    INFO("Terms with 2 operators");
    L.printTerms(2);

    LatticePresets::addSS(&L, "K", chain_site_names[0], J);
    INFO("Terms with 4 operators");
    L.printTerms(4);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare(true);
    print_section("Indices");
    IndexInfo.printIndices();

    print_section("Matrix element storage");
    IndexHamiltonian Storage(&L,IndexInfo);
    Storage.prepare();
    print_section("Terms");
    INFO(Storage);

    Symmetrizer Symm(IndexInfo, Storage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();
    
    auto Q=Symm.getQuantumNumbers();
    Q.set(0,2);
    Q.set(1,0.0);
    INFO(Q);
    auto B=S.getBlockNumber(Q);
    INFO("Looking block " << B << " with quantum numbers " << Q);
    HamiltonianPart Hpart(IndexInfo, Storage, S, B);
    Hpart.prepare();
    INFO_NONEWLINE("Diagonalizing...");
    Hpart.diagonalize();
    INFO("done.");
    size_t IndexSize = IndexInfo.getIndexSize();
    
    Operator Splus;
    Operator Sminus;
    std::vector<ParticleIndex> up_indices, down_indices;
    for (auto site_pairs : L.getSiteMap()) {  
        auto site_name = site_pairs.first;
        ParticleIndex up = IndexInfo.getIndex(site_name, 0, 0);
        ParticleIndex down = IndexInfo.getIndex(site_name, 0, 1);
        up_indices.push_back(up);
        down_indices.push_back(down);
        Splus += OperatorPresets::Cdag(up)*OperatorPresets::C(down);
        Sminus  += OperatorPresets::Cdag(down)*OperatorPresets::C(up);
        }
    Operator Sz = OperatorPresets::Sz(up_indices,down_indices);
    Operator SzSz = Sz*Sz;
    Operator SplusSminus = Splus*Sminus;
    Operator SminusSplus = Sminus*Splus;
    Operator S2 = SzSz + (SplusSminus + SminusSplus)*0.5;
    Operator S2_2 = SzSz + SplusSminus - Sz;
    Operator S2_3 = SzSz + SminusSplus + Sz;

    Operator Sc = Splus.getCommutator(Sz);
    DEBUG((Sc == Splus *(-1.0)));
    Operator Sc2 = Sminus.getCommutator(Sz);
    DEBUG((Sc2 == Sminus ));
    DEBUG((S2 == S2_2));
    DEBUG((S2 == S2_3));
    DEBUG(S2.commutes(Sz));
    DEBUG(S2.commutes(Splus));

    auto blockstates = S.getFockStates(B);
    
    for (size_t state_index = 0; state_index < 1; ++state_index) {
        auto State = Hpart.getEigenState(state_index);

        for (size_t i=0; i<State.size(); ++i) { if (!is_equal(State(i),0.0,1e-3)) INFO_NONEWLINE(State(i) << "*|" << blockstates[i] << "> + "); }; INFO("");
        //RealType s2val = S2.getMatrixElement(State,State,blockstates);
        //RealType szszval = SzSz.getMatrixElement(State,State,blockstates); 
        MelemType splussminusval = SplusSminus.getMatrixElement(State,State,blockstates);

        //INFO("<S^2> = " << s2val);
        //INFO("<SzSz> = " << szszval);
        INFO("<S+S-> = " << splussminusval);
    };

/*
    Hamiltonian H(IndexInfo, Storage, S);
    H.prepare();
    H.diagonalize();
    INFO("The value of ground energy is " << H.getGroundEnergy());

    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(H.getEigenValue(i)); 
 
    //srand (time(NULL));

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

    std::ofstream F("Gw00.dat");
    for(int n = 0; n<wn_max; ++n) {
        ComplexType GF_val = GF(n);
        F << real(GF_val) << "  " << imag(GF_val) << std::endl;
        }
    F.close();
*/
}
