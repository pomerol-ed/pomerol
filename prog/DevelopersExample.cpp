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
/*
#include "LatticeAnalysis.h"
#include "Term.h"
#include "IndexClassification.h"
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

void printFramed (const std::string& str)
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

  JSONLattice L;
  //L.readin(opt.LatticeFile);
  Lattice::TermStorage s1; 

  Lattice::Term T1(2);
  T1.Order[0]=false;
  T1.Order[1]=true;
  T1.Sites[0]="A";
  T1.Sites[1]="B";
  T1.Spins[0]=0;
  T1.Spins[1]=0;
  T1.Orbitals[0]=0;
  T1.Orbitals[1]=0;
  T1.Value=1.0;

  cout << T1 << endl;
  Lattice::Term *T=&T1;
  s1.addTerm(T);
  cout << T->getOrder() << endl;

  Lattice::Term T2(4);
  T2.Order[0]=false;
  T2.Order[1]=true;
  T2.Order[2]=false;
  T2.Order[3]=true;
  T2.Sites[0]="A";
  T2.Sites[1]="B";
  T2.Sites[2]="C";
  T2.Sites[3]="D";

  T=&T2;
  cout << T2 << endl;
  cout << T->getOrder() << endl;
  s1.addTerm(T);

  s1.getTermList(2);

  Lattice::Presets::addSSite((Lattice*) &L,  std::string("A"), 1.0, 0.5);  
  Lattice::Presets::addSSite((Lattice*) &L,  std::string("B"), 2.0, 0.5, 3, 3);  

  return 0;
};


