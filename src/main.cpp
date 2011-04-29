#include "Misc.h"
#include "output.h"
#include "HDF5Storage.h"
#include "LatticeAnalysis.h"
#include "Term.h"
#include "IndexClassification.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "FieldOperator.h"
#include "GFContainer.h"
#include "TwoParticleGFContainer.h"
#include "Vertex4.h"

#include "OptionParser.h"

#include <fstream>
std::string input = "system.ini";

LatticeAnalysis Lattice;
IndexClassification IndexInfo(Lattice);
StatesClassification S(IndexInfo); 
output_handle OUT;
Hamiltonian H(IndexInfo,S,OUT,input);

std::ostream &OUTPUT_STREAM=std::cout;


/* ======================================================================== */
// To be removed

void printFramed (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

enum AmpStyle{UnAmputated, Amputated};

inline RealType chop(RealType &i){ return (std::fabs(i)<1e-10)?0.0:i; }

void saveChi(const char *fname, TwoParticleGFContainer &Chi, int size_wg) 
{
  std::cout << "Dumping Chi4..." << std::flush;
#if defined(HRD)
  std::ofstream chi_str(fname,std::ios::out);
  chi_str.setf(std::ios::scientific, std::ios::floatfield);
  chi_str.setf(std::ios::showpoint);
  chi_str.precision(8);
  chi_str<<"Re              Im                       z1 z2          w1' w1 w2' w2           n1' n1 n2' n2\n";
#else
  std::ofstream chi_str(fname,std::ios::out | std::ios::binary);
#endif
  int n_zone=2;
  int n_part=IndexInfo.getIndexSize()/2.;
  RealType acc=1e-8;
  for (int z1=0; z1<n_zone; z1++)
    for (int z2=n_zone-1; z2>=0; z2--)
      for(int n1=0; n1<n_part; n1++)
        for(int n1_=0; n1_<n_part; n1_++)
          for(int n2=0; n2<n_part; n2++)
            for(int n2_=0; n2_<n_part; n2_++)
              for(int w1=-size_wg; w1<size_wg; w1++)
                for(int w1_=-size_wg; w1_<size_wg; w1_++)
                  for(int w2=-size_wg; w2<size_wg; w2++){

                    int w2_=w1+w2-w1_;
                    if (w2_>=-size_wg && w2_<size_wg){
                      TwoParticleGFContainer::IndexCombination *comb1;
                      comb1 = new TwoParticleGFContainer::IndexCombination(n1+n_part*z1,n2+n_part*z2,n1_+z1*n_part,n2_+z2*n_part);
                      std::complex<double> *z=new ComplexType(Chi(*comb1,w1,w2,w1_));
                      delete comb1;
                      if(abs(*z)>acc){
#if defined(HRD)
                        chi_str << chop(real(*z)) <<"  "<< chop(imag(*z)) << "           "
                          << z1 <<" "<< z2 << "           " << w1 << "  " << w1_ << " " << w2 << "  " << w2_ 
                          << "            "<<n1<<"  "<<n1_<<" "<<n2<<"  "<<n2_<< "            "
                          << std::endl << std::flush;
#else
                        chi_str.write(reinterpret_cast<char *>(z),sizeof(std::complex<double>));

                        chi_str.write(reinterpret_cast<char *>(&z1),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&z2),sizeof(int));

                        chi_str.write(reinterpret_cast<char *>(&w1),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&w1_),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&w2),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&w2_),sizeof(int));

                        chi_str.write(reinterpret_cast<char *>(&n1),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&n1_),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&n2),sizeof(int));
                        chi_str.write(reinterpret_cast<char *>(&n2_),sizeof(int));

#endif
                      };
                      delete z;
                    }
                  }
  chi_str<<"0 0"<<std::endl;
  std::cout << "Finished." << std::endl;
  return;
}			    

void saveGamma(const char *fname, Vertex4 &Vertex, std::vector<TwoParticleGFContainer::IndexCombination*>& Combinations, int size_wg, unsigned short style = Amputated) 
{
  std::cout << "Dumping Gamma4..." << std::flush;
#if defined(HRD)
  std::ofstream gamma_str(fname,std::ios::out);
  gamma_str.setf(std::ios::scientific, std::ios::floatfield);
  gamma_str.setf(std::ios::showpoint);
  gamma_str.precision(8);
  gamma_str<<"Re              Im                       z1 z2          w1' w1 w2' w2           n1' n1 n2' n2\n";
#else
  std::ofstream gamma_str(fname,std::ios::out | std::ios::binary);
#endif
  int n_zone=2;
  int n_part=IndexInfo.getIndexSize()/2.;
  RealType acc=1e-10;
  int z1,z2,n1,n1_,n2,n2_;
  for (std::vector<TwoParticleGFContainer::IndexCombination*>::const_iterator comb1=Combinations.begin(); comb1!=Combinations.end(); ++comb1)
              for(int w1=-size_wg; w1<size_wg; w1++)
                for(int w1_=-size_wg; w1_<size_wg; w1_++)
                  for(int w2=-size_wg; w2<size_wg; w2++){

                    int w2_=w1+w2-w1_;
                    if (w2_>=-size_wg && w2_<size_wg){
                      ComplexType *z;
                      
                      //comb1 = TwoParticleGFContainer::IndexCombination(n1+n_part*z1,n2+n_part*z2,n1_+z1*n_part,n2_+z2*n_part);
                      n1 =  (*comb1)->Indices[0]%n_part;
                      n2 =  (*comb1)->Indices[1]%n_part;
                      n1_ = (*comb1)->Indices[2]%n_part;
                      n2_ = (*comb1)->Indices[3]%n_part;
                      z1 = (*comb1)->Indices[0]/n_part;
                      z2 = (*comb1)->Indices[1]/n_part;
                      z=new ComplexType(Vertex(**comb1,w1,w2,w1_));
                      if(abs(*z)>acc){
#if defined(HRD)
                        gamma_str << chop(real(*z)) <<"  "<< chop(imag(*z)) << "           "
                          << z1 <<" "<< z2 << "           " << w1 << "  " << w1_ << " " << w2 << "  " << w2_ 
                          << "            "<<n1<<"  "<<n1_<<" "<<n2<<"  "<<n2_<< "            "
                          << std::endl << std::flush;
#else
                        gamma_str.write(reinterpret_cast<char *>(z),sizeof(std::complex<double>));

                        gamma_str.write(reinterpret_cast<char *>(&z1),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&z2),sizeof(int));

                        gamma_str.write(reinterpret_cast<char *>(&w1),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&w1_),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&w2),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&w2_),sizeof(int));

                        gamma_str.write(reinterpret_cast<char *>(&n1),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&n1_),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&n2),sizeof(int));
                        gamma_str.write(reinterpret_cast<char *>(&n2_),sizeof(int));

#endif
                      };
                     delete z;
                    }
                  }
  gamma_str<<"0 0"<<std::endl;
  std::cout << "Finished." << std::endl;
  return;
}			    

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

  printFramed("Lattice Info");
  Lattice.readin(opt.LatticeFile);
  std::cout << Lattice.printSitesList().str() << std::flush;
  IndexInfo.prepare();
  printFramed("System Info");
  IndexInfo.printIndexInfoList();
  printFramed("Hopping Matrix");
  IndexInfo.printHoppingMatrix();
  printFramed("Terms check");
  IndexInfo.printTerms();
  printFramed("Equivalent Permutations");
  IndexInfo.printEquivalentPermutations();
  OUT = output_handle("output");

  S.iniStatesClassification();

  //end of determination    

  printFramed("System is determined");
  printFramed("Process of creation and diagonalization all parts of Hamiltonian has started");
  
  //begining of creation all part of Hammiltonian

  H.enter();
  H.dump();
  H.diagonalize();
  RealType beta = opt.beta;
  H.dump();

  num_cout << std::endl << "The value of ground energy is " << H.getGroundEnergy() << std::endl;

 //   GFContainer G(S,H,rho,IndexInfo,Operators);
 //   G.prepare();
 //   G.compute();
 //   cout << std::setprecision(9) << G(0) << endl;
 //   cout << std::setprecision(9) << G(1) << endl;

  DensityMatrix rho(S,H,beta);
  //    DensityMatrix.reduce();
  rho.prepare();
  rho.compute();
  num_cout << "<H> = " << rho.getAverageEnergy() << std::endl;

  try{
    HDF5Storage dmp1("test1.h5");
    dmp1.save(rho);
  
    DensityMatrix rho_loaded(S,H,beta);
    dmp1.load(rho_loaded);
  
    HDF5Storage dmp2("test2.h5");
    dmp2.save(rho_loaded);
  } catch(H5::Exception& ex){
      ex.printError();
  }


  /*   for (QuantumState i=0; i < S.N_st(); ++i) 
       cout << std::setw(20) << "E:" << H.eigenval(i) << "\t E-E0 " << H.eigenval(i) - rho.getAverageEnergy() << "\t weight: " << rho(i) << "  " << exp(-beta*(H.eigenval(i) - H.getGroundEnergy()))/1.000 << endl; 
       */
  std::ofstream out;
  out.open("output/Stat.En.dat");
  out << iomanip_prefs << rho.getAverageEnergy() << std::endl;
  out.close();

  out.open("output/Stat.NN.dat");
  out << iomanip_prefs << rho.getAverageDoubleOccupancy(0,IndexInfo.getIndexSize()/2.) << std::endl;
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
    G.dumpToPlainText(2*wn);

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
    saveGamma("Gamma4.dat",Gamma4,v1,wn,Amputated);

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
  return 0;
};


