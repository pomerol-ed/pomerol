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
#include "Operator.h"
#include "IndexHamiltonian.h"
#include "Symmetrizer.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "FieldOperatorPart.h"
#include "FieldOperator.h"
#include "FieldOperatorContainer.h"
#include "DensityMatrixPart.h"
#include "DensityMatrix.h"
#include "GFContainer.h"
#include "TwoParticleGFContainer.h"

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

  print_section("Lattice");
  INFO("Sites");
  L->printSites();
  INFO("Terms with 2 operators");
  L->printTerms(2);
  INFO("Terms with 4 operators");
  L->printTerms(4);

  IndexClassification IndexInfo(L->getSiteMap());
  IndexClassification *II = &IndexInfo;
  IndexInfo.prepare();
  print_section("Indices");
  IndexInfo.printIndices();

  print_section("Matrix element storage");
  IndexHamiltonian Storage(L,IndexInfo);
  Storage.prepare();
  INFO("Terms with 2 operators");
  Storage.printTerms(2);
  INFO("Terms with 4 operators");
//  Storage.printTerms(4);

  //DEBUG("Check - all terms");
  //Storage.printAllTerms();
  Symmetrizer Symm(IndexInfo, Storage);
  Symm.compute();
  StatesClassification S(IndexInfo,Symm);

  S.compute();

  Hamiltonian H(IndexInfo, Storage, S);
  H.prepare();
  H.diagonalize();
  INFO("The value of ground energy is " << H.getGroundEnergy());
  
  FieldOperatorContainer Operators(IndexInfo, S, H);
  Operators.prepare();

  RealType beta=opt.beta;
  DensityMatrix rho(S,H,beta);

  rho.prepare();
  rho.compute();
  INFO("<H> = " << rho.getAverageEnergy() << std::endl)
  GFContainer G(IndexInfo,S,H,rho,Operators);
  long wn = opt.NumberOfMatsubaras;


  if (1==1) {

      print_section("Green's function calculation");
      std::set<IndexCombination2> v2;
      v2.insert(IndexCombination2(0,0));
      v2.insert(IndexCombination2(IndexInfo.getIndexSize()/2,IndexInfo.getIndexSize()/2));
      G.prepareAll(v2);
      G.computeAll();
       }

  if (1==1) {   
      print_section("Two Particle Green's function calculation");
      std::set<IndexCombination4> v1;
      v1.insert(IndexCombination4(0,IndexInfo.getIndexSize()/2,0,IndexInfo.getIndexSize()/2));
      v1.insert(IndexCombination4(0,0,0,0));
      TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
      Chi4.prepareAll(v1);
      Chi4.computeAll(wn);
      }

  return 0;
};


