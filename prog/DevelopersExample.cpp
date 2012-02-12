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
#include "LatticeReader.h"
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
    
  LatticeReader L;
  L.readinFromJSON(opt.LatticeFile);
  Json::Value A=L.getDictionary();
  cout << A["General"]["mu"].asDouble() << endl;
  cout << A["General"] << endl;
  exit(0);

  LatticeAnalysis Lattice;
  IndexClassification IndexInfo(Lattice);
  StatesClassification S(IndexInfo); 
  Hamiltonian H(IndexInfo,S);

  printFramed("Lattice Info");
  Lattice.readin(opt.LatticeFile);
  std::cout << Lattice.printSitesList().str() << std::flush;
  IndexInfo.prepare();
  printFramed("System Info");
  IndexInfo.printIndexList();
  printFramed("Hopping Matrix");
  IndexInfo.printHoppingMatrix();
  printFramed("Terms check");
  IndexInfo.printTerms();
  printFramed("Equivalent Permutations");
  IndexInfo.printEquivalentPermutations();

  S.compute();

  //end of determination    

  printFramed("System is determined");
  printFramed("Process of creation and diagonalization all parts of Hamiltonian has started");
  
  exit(0);
  #warning Something is wrong with HDF5Storage. segfaults detected in this rev
  //  HDF5Storage storage("test.h5"); 
  
  //begining of creation all part of Hammiltonian

  H.prepare();
  H.compute();
  RealType beta = opt.beta;
  //H.dump();
  //storage.save(H);

   INFO("The value of ground energy is " << H.getGroundEnergy());

 //   GFContainer G(S,H,rho,IndexInfo,Operators);
 //   G.prepare();
 //   G.compute();
 //   cout << std::setprecision(9) << G(0) << endl;
 //   cout << std::setprecision(9) << G(1) << endl;

  DensityMatrix rho(S,H,beta);
  //    DensityMatrix.reduce();
  rho.prepare();
  rho.compute();
  INFO("<H> = " << rho.getAverageEnergy() << std::endl)

  /*   for (QuantumState i=0; i < S.N_st(); ++i) 
       cout << std::setw(20) << "E:" << H.eigenval(i) << "\t E-E0 " << H.eigenval(i) - rho.getAverageEnergy() << "\t weight: " << rho(i) << "  " << exp(-beta*(H.eigenval(i) - H.getGroundEnergy()))/1.000 << endl; 
       */
  std::ofstream out;
  out.open("output/Stat.En.dat");
  out << rho.getAverageEnergy() << std::endl;
  out.close();

  out.open("output/Stat.NN.dat");
  out << rho.getAverageDoubleOccupancy(0,IndexInfo.getIndexSize()/2.) << std::endl;
  out.close();

  //finishing of creation
  std::cout << std::endl;
  std::cout << "All parts are created!" << std::endl;
  std::cout << std::endl;

  FieldOperatorContainer Operators(S,H,IndexInfo);
  GFContainer G(S,H,rho,IndexInfo,Operators);
  long wn = opt.NumberOfMatsubaras;

  if (1==1){
    printFramed("Two Particle Green's function calculation");

    std::vector<GFContainer::IndexCombination*> v2;
    v2.push_back(new GFContainer::IndexCombination(0,0));
    v2.push_back(new GFContainer::IndexCombination(IndexInfo.getIndexSize()/2,IndexInfo.getIndexSize()/2));
    G.readInitialIndices(v2);
    G.prepare();
    G.compute();
    G.dumpToPlainText(8*wn);

    std::vector<TwoParticleGFContainer::IndexCombination*> v1;
    v1.push_back(new TwoParticleGFContainer::IndexCombination(0,IndexInfo.getIndexSize()/2,0,IndexInfo.getIndexSize()/2));
    v1.push_back(new TwoParticleGFContainer::IndexCombination(0,0,0,0));
    TwoParticleGFContainer Chi4(S,H,rho,IndexInfo,Operators);
    Chi4.readInitialIndices(v1);
    Chi4.prepare();
    Chi4.compute(wn);
    //saveChi("Chi4.dat",Chi4,wn);
    
    Vertex4 Gamma4(IndexInfo,Chi4,G);
    Gamma4.prepareUnAmputated();
    Gamma4.computeUnAmputated();
    Gamma4.prepareAmputated(v1);
    Gamma4.computeAmputated();

/*
    for (unsigned short i=0;i<v1.size();i++){
      std::cout << "Chi4" << *v1[i] << std::endl;
      std::cout << Chi4(*v1[i],3,2,0) << std::endl;
      std::cout << Chi4(*v1[i],2,5,2) << std::endl;
      std::cout << Chi4(*v1[i],5,2,5) << std::endl;
      std::cout << Chi4(*v1[i],5,2,2) << std::endl;
      std::cout << Chi4(*v1[i],1,7,1) << std::endl;
      std::cout << Chi4(*v1[i],2,-2,4) << std::endl;
      std::cout << Chi4(*v1[i],29,-29,29) << std::endl << std::endl;
    };

    for (unsigned short i=0;i<v1.size();i++){
      std::cout << "Gamma4" << *v1[i] << std::endl;
      std::cout << Gamma4(*v1[i],3,2,0) << std::endl;
      std::cout << Gamma4(*v1[i],2,5,2) << std::endl;
      std::cout << Gamma4(*v1[i],5,2,5) << std::endl;
      std::cout << Gamma4(*v1[i],5,2,2) << std::endl;
      std::cout << Gamma4(*v1[i],1,7,1) << std::endl;
      std::cout << Gamma4(*v1[i],2,-2,4) << std::endl;
      std::cout << Gamma4(*v1[i],29,-29,29) << std::endl << std::endl;
 
    };
*/
    //   comb1 = new TwoParticleGFContainer::IndexCombination(0,2,0,1);
    //   cout << Chi4.vanishes(*comb1) << endl;
    //DEBUG(Chi4.getNumNonResonantTerms() << " non-resonant terms");
    //DEBUG(Chi4.getNumResonantTerms() << " resonant terms");    

  };
 // DEBUG HDF5 save/load
  //storage.close();
  //HDF5Storage storage_load("test.h5");
  //Hamiltonian H2(IndexInfo,S,OUT,input);
  //storage_load.load(H2);
  //HDF5Storage storage2("test2.h5");
  //storage2.save(H2);
  //exit(0);
  //RowMajorMatrixType SP1(2,3);
  //SP1.startVec(0);
  //SP1.insertBack(0,0) = 1;
  //SP1.insertBack(0,1) = 2;
  //SP1.insertBack(0,2) = 3;
  //SP1.startVec(1);
  //SP1.insertBack(1,0) = 4;
  //SP1.insertBack(1,1) = 5;
  //SP1.insertBack(1,2) = 6;
  //SP1.finalize();
  
  //DEBUG("SP1" << SP1)
  //HDF5Storage::saveRowMajorMatrix(&storage,"SP",SP1);
  //RowMajorMatrixType SP2;
  //HDF5Storage::loadRowMajorMatrix(&storage,"SP",SP2);
  //DEBUG("SP2" << SP2)
  //exit(0);


  return 0;
};


