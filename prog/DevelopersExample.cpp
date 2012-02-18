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

  Lattice::Presets::addSSite((Lattice*) &L,  std::string("A"), 1.0, 0.5);  
  Lattice::Presets::addSSite((Lattice*) &L,  std::string("B"), 2.0, 0.5);  
  L.addTerm(Lattice::Term::Presets::Hopping(std::string("A"), std::string("B"), 1.0, 0, up));
  L.addTerm(Lattice::Term::Presets::Hopping(std::string("A"), std::string("B"), 1.0, 0, down));
  Lattice::Presets::addPSite((Lattice*) &L,  std::string("C"), 4, 1, 5, 2, 2);  

  INFO("Terms with 2 operators");
  L.printTerms(2);
  INFO("Terms with 4 operators");
  L.printTerms(4);

  return 0;
};


