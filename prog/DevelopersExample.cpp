//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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


#include "Misc.h"
#include "HDF5Storage.h"
#include "Lattice.h"
#include "LatticePresets.h"
#include "IndexClassification.h"
#include "IndexHamiltonian.h"
#include "Symmetrizer.h"
/*
#include "LatticeAnalysis.h"
#include "Term.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "FieldOperator.h"
#include "GFContainer.h"
#include "TwoParticleGFContainer.h"
#include "Vertex4.h"
#include "Logger.h"
*/
#include "OptionParser.h"

#include <fstream>

using namespace Pomerol;
using std::string; using std::cout; using std::endl;

/* ======================================================================== */
// To be removed

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

enum AmpStyle{UnAmputated, Amputated};

/* ======================================================================== */

int main(int argc, char *argv[])
{
  pomerolOptionParser opt;
   try {
		opt.parse(&argv[1], argc-1); // Skip argv[0].

        std::cout << "pomerolDiag. Parameters " << std::endl;
		std::cout << "Lattice File         : " << opt.LatticeFile << std::endl;
		std::cout << "Number Of Matsubaras : " << opt.NumberOfMatsubaras << std::endl;
		std::cout << "beta:                : " << opt.beta << std::endl;
	} catch (const optparse::unrecognized_option& e) {
		std::cout << "unrecognized option: " << e.what() << std::endl;
		return 1;
	} catch (const optparse::invalid_value& e) {
		std::cout << "invalid value: " << e.what() << std::endl;
		return 1;
	}

  Log.setDebugging(true);

  JSONLattice JL;
  Lattice *L=&JL;
  JL.readin(opt.LatticeFile);

  /* std::vector<ParticleIndex> a1(3);
  a1[0]=1;
  a1[1]=2;
  a1[2]=0;
  Symmetrizer::IndexPermutation AA2(a1);
  DynamicIndexCombination AA(3);
  AA[0]=3;
  AA[1]=3;
  AA[2]=3;
  DynamicIndexCombination BB(3);
  BB[0]=3;
  BB[1]=3;
  BB[2]=3;
  DEBUG(AA2<"< " << BB << " = " << (AA<BB));
  DEBUG(AA<<"==" << BB << " = " << (AA==BB));
  DEBUG(AA<<"!=" << BB << " = " << (AA!=BB));
  exit(0);

  L->addSite(new Lattice::Site("A",1,2));
  L->addSite(new Lattice::Site("B",2,2));
  L->addSite(new Lattice::Site("C",1,2));
  Lattice::Presets::addHopping(L,"A","B",1.37,0,1);
  Lattice::Presets::addHopping(L,"B","C",1.0,1,0);
  Lattice::Presets::addCoulombP(L, "B",4, 0.5, 0.0);
  Lattice::Presets::addHopping(L,"A","C",0.612);
  */
  print_section("Lattice");
  INFO("Sites");
  L->printSites();
  INFO("Terms with 2 operators");
  L->printTerms(2);
  INFO("Terms with 4 operators");
  L->printTerms(4);

    //DEBUG(L->getTermStorage().getMaxTermOrder());
  IndexClassification Indices(L->getSiteMap());
  IndexClassification *II = &Indices;
  Indices.prepare();
  print_section("Indices");
  Indices.printIndices();
  //DEBUG(Indices.getIndex("B",0,1));
  //DEBUG(Indices.checkIndex(Indices.getIndex("B",0,3)));

  print_section("Matrix element storage");
  IndexHamiltonian Storage(L,Indices);
  Storage.prepare();
  INFO("Terms with 2 operators");
  Storage.printTerms(2);
  INFO("Terms with 4 operators");
  Storage.printTerms(4);

  Symmetrizer S(Indices, Storage);
  return 0;
};


