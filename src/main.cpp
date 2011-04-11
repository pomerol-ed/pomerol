#define EIGEN2_SUPPORT

#ifdef DMTruncate
#define DENSITY_MATRIX_TRUNCATION_TOLERANCE 1e-3
#define TERM_MATRIX_ELEMENT_TOLERANCE 1e-3
#define TERM_RESONANCE_TOLERANCE 1e-16
#endif

#include "config.h"
#include "output.h"
#include "LatticeAnalysis.h"
#include "Term.h"
#include "BitClassification.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "FieldOperatorPart.h"
#include "FieldOperator.h"
#include "GreensFunction.h"
#include "GFContainer.h"
#include "TwoParticleGF.h"
#include "TwoParticleGFContainer.h"
#include "Vertex4.h"

#include <iostream>
#include <fstream>

#include "iniconfig.h"

string input = "system.ini";

LatticeAnalysis Lattice;
BitClassification IndexInfo(Lattice);
StatesClassification S(IndexInfo); 
output_handle OUT;
Hamiltonian H(IndexInfo,S,OUT,input);

ostream &OUTPUT_STREAM=std::cout;

IniConfig* pIni;

/* ======================================================================== */
// To be removed

RealType chop(RealType &i){ return (std::fabs(i)<1e-5)?0.0:i; }

void saveChi(const char *fname, TwoParticleGFContainer &Chi, int size_wg) 
{
  std::cout << "Dumping Chi4..." << std::flush;
  std::ofstream chi_str(fname,std::ios::out);
  chi_str.setf(std::ios::fixed, std::ios::floatfield);
  chi_str.setf(std::ios::showpoint);
  chi_str.precision(8);
  chi_str<<"Re         Im               z1 z2          w1' w1 w2' w2        n1' n1 n2' n2\n";
  int n_zone=2;
  int n_part=S.N_b()/2.;
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
                        //std::cout<<z1<<" "<<z2<<" "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<" "<<w1<<" "<<w1_<<" "<<w2<<" "<<w2_<<" bb "<<std::endl;
                        std::complex<double> z=Chi(*comb1,w1,w2,w1_);
                        //std::cout<<z1<<" "<<z2<<" "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<" "<<w1<<" "<<w1_<<" "<<w2<<" "<<w2_<<" "<<z<<std::endl;
                        if(abs(z)>acc){
                            //int w1_=w1+W, w2=w2_+W;
                            chi_str << chop(real(z)) <<"  "<< chop(imag(z)) << "           "
                            << z1 <<" "<< z2 << "           " << w1 << "  " << w1_ << " " << w2 << "  " << w2_ 
                            << "            "<<n1<<"  "<<n1_<<" "<<n2<<"  "<<n2_<< "            "
                            << endl << flush;
                        }
                    }
                  }
  chi_str<<"0 0"<<std::endl;
  std::cout << "Finished." << std::endl;
  return;
}			    
/* ======================================================================== */

int main()
{   
        
    cout << "=======================" << endl;
    cout << "Lattice Info" << endl;
    cout << "=======================" << endl;
    Lattice.readin();
    cout << Lattice.printSitesList().str() << flush;
    IndexInfo.prepare();
    cout << "=======================" << endl;
    cout << "System Info" << endl;
    cout << "=======================" << endl;
    IndexInfo.printBitInfoList();
    cout << "=======================" << endl;
    cout << "Hopping Matrix" << endl;
    cout << "=======================" << endl;
    IndexInfo.printHoppingMatrix();
    cout << "=======================" << endl;
    cout << "Terms check" << endl;
    cout << "=======================" << endl;
    IndexInfo.printTerms();
    //determination of system
    
    pIni = new IniConfig(input);
    
    OUT = output_handle((*pIni)["output:path"]);

    S.iniStatesClassification();

    //end of determination    
    
    cout << "=======================" << endl;
    cout << "System is determined" << endl;
    cout << "=======================" << endl;
    cout << "=======================================" << endl;
    cout << "Process of creation and diagonalization" << endl;
    cout << "all parts of  Hamiltonian  has  started" << endl;
    cout << endl;

    //begining of creation all part of Hammiltonian

    H.enter();
    H.dump();
    H.diagonalize();
    RealType beta = (*pIni)["Green Function:beta"];
    RealType ProbabilityCutoff = RealType((*pIni)["system:ProbabilityCutoff"]);
    if (ProbabilityCutoff>0 && ProbabilityCutoff < 1) H.reduce(-log(ProbabilityCutoff)/beta);
    H.dump();

    num_cout << endl << "The value of ground energy is " << H.getGroundEnergy() << endl;

   
    DensityMatrix rho(S,H,beta);
//    DensityMatrix.reduce();
    rho.prepare();
    rho.compute();
    num_cout << "<H> = " << rho.getAverageEnergy() << endl;

    
 /*   for (QuantumState i=0; i < S.N_st(); ++i) 
        cout << std::setw(20) << "E:" << H.eigenval(i) << "\t E-E0 " << H.eigenval(i) - rho.getAverageEnergy() << "\t weight: " << rho(i) << "  " << exp(-beta*(H.eigenval(i) - H.getGroundEnergy()))/1.000 << endl; 
        */
    std::ofstream out;
    out.open("output/Stat.En.dat");
    out << iomanip_prefs << rho.getAverageEnergy() << endl;
    out.close();
    
      


    //finishing of creation
    cout << endl;
    cout << "All parts are created!" << endl;
    cout << endl;

    FieldOperatorContainer Operators(S,H,IndexInfo);
    GFContainer G(S,H,rho,IndexInfo,Operators);
    G.prepare();
    G.compute();
    cout << std::setprecision(9) << G(0) << endl;
    cout << std::setprecision(9) << G(1) << endl;

    if ((*pIni)["System:calculate_2PGF"]){
        cout << endl;
        cout << "==========================================" << endl;
        cout << "Two Particle Green's function calculation" << endl;
        cout << "==========================================" << endl;


        TwoParticleGFContainer::IndexCombination *comb1;
        std::vector<TwoParticleGFContainer::IndexCombination*> v1;

        comb1 = new TwoParticleGFContainer::IndexCombination(0,0,0,0);
        v1.push_back(comb1);
        comb1 = new TwoParticleGFContainer::IndexCombination(0,1,0,1);
        v1.push_back(comb1);
//
        int wn = (int) (*pIni)["Green Function:points"];
        TwoParticleGFContainer Chi4(S,H,rho,IndexInfo,Operators);
        Chi4.readNonTrivialIndices(v1);
        Chi4.prepare();
        Chi4.compute(wn);

        saveChi("Chi4.dat",Chi4,wn);

        for (unsigned short i=0;i<v1.size();i++){
            cout << std::setprecision(9) << Chi4(*v1[i],3,2,0) << endl;
            cout << Chi4(*v1[i],2,5,2) << endl;
            cout << Chi4(*v1[i],5,2,2) << endl;
            cout << Chi4(*v1[i],1,7,1) << endl;
            cout << Chi4(*v1[i],2,-2,4) << endl;
            cout << Chi4(*v1[i],29,-29,29) << endl << endl;
        cout << *v1[i] << " : " << (bool) Chi4.vanishes(*v1[i]) << endl;
        };
        
     //   comb1 = new TwoParticleGFContainer::IndexCombination(0,2,0,1);
     //   cout << Chi4.vanishes(*comb1) << endl;
        //DEBUG(Chi4.getNumNonResonantTerms() << " non-resonant terms");
        //DEBUG(Chi4.getNumResonantTerms() << " resonant terms");    
    };
};


