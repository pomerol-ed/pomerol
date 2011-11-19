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
#include <H5File.h>
#include <H5Group.h>
std::string input = "system.ini";

LatticeAnalysis Lattice;
IndexClassification IndexInfo(Lattice);
StatesClassification S(IndexInfo); 
Hamiltonian H(IndexInfo,S);

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

// An UGLY workaround! We have to reorganize SingleIndex class and
// incorporate this method into it.
std::string indexToString(SingleIndex* index)
{
    std::stringstream strs;

    strs << "(";
    strs << index->site;
    strs << "," << (index->spin == 0 ? "down" : "up");
    strs << "," << index->type;
    switch(index->type){
	case 1:
	    strs << "," << reinterpret_cast<pSingleIndex*>(index)->orbital;
	    break;
	//case 2:
	//break;
	default:
	    strs << ",0";
	    break;
    }
    strs << ")";

    return strs.str();
}

void addServiceInfo(H5::H5File& file)
{
    // state_type
    H5::EnumType state_type(H5::PredType::NATIVE_SHORT);
    int state_type_values[2] = {0,1};
    state_type.insert("CREATE",&state_type_values[0]);
    state_type.insert("PLACEHOLDER",&state_type_values[1]);
    state_type.commit(file,"state_type");

    // log_type
    H5::CompType log_type(2*sizeof(char*));
    H5::StrType VarCStr(H5::PredType::C_S1);
    VarCStr.setSize(H5T_VARIABLE);
    log_type.insertMember("time",0,VarCStr);
    log_type.insertMember("log",sizeof(char*),VarCStr);
    log_type.commit(file,"log_type");

    // revisions
    H5::Group revisions_group = file.createGroup("revisions");
    H5::Attribute revisions_last_attr = revisions_group.createAttribute("last",H5::PredType::NATIVE_INT,H5::DataSpace());
    int zero = 0;
    revisions_last_attr.write(H5::PredType::NATIVE_INT,&zero);
}

void complexify(H5::DataSet& set)
{
    // Add "__complex__" attribute
    H5::Attribute complex_attr = set.createAttribute("__complex__",H5::PredType::NATIVE_CHAR,H5::DataSpace());
    char _true = 1;
    complex_attr.write(H5::PredType::NATIVE_CHAR,&_true);
}

void writeGForHartmut(GFContainer& G, long NumberOfMatsubaras)
{
    std::vector<SingleIndex*> Indices1 = IndexInfo.getSingleIndexList();
    std::vector<SingleIndex*> Indices2 = IndexInfo.getSingleIndexList();

    hsize_t Dims[2] = {NumberOfMatsubaras,2};
    H5::DataSpace DS(2,Dims);
    
    for(std::vector<SingleIndex*>::iterator pIndex1 = Indices1.begin();pIndex1 != Indices1.end();pIndex1++)
    for(std::vector<SingleIndex*>::iterator pIndex2 = Indices2.begin();pIndex2 != Indices2.end();pIndex2++){
	
	std::string fileName = "g_" + indexToString(*pIndex1) + indexToString(*pIndex2) + ".h5";
	H5::H5File hdf5file(fileName,H5F_ACC_EXCL); // Do not overwrite existing files!
	
	addServiceInfo(hdf5file);
	H5::DataSet DataSet = hdf5file.createDataSet("/data",H5::PredType::NATIVE_DOUBLE,DS);
	double data[2*NumberOfMatsubaras];
	for(long wn = 0; wn < NumberOfMatsubaras; ++wn){
	    std::complex<double> value;
	    value = G((*pIndex1)->bitNumber,(*pIndex2)->bitNumber,wn);
	  
	    data[2*wn] = real(value);
	    data[2*wn+1] = imag(value);
	}
	DataSet.write(data,H5::PredType::NATIVE_DOUBLE);
	complexify(DataSet);
		
	hdf5file.close();
    }
}

/*inline RealType chop(RealType &i){ return (std::fabs(i)<1e-10)?0.0:i; }

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
                      z=new ComplexType(Vertex(**comb1,w1,w2,w1_)*(-1.));
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

    std::cout << "Finished." << std::endl;
    return;
}			    
*/
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
   
    //begining of creation all part of Hammiltonian

    H.prepare();
    H.compute();
    RealType beta = opt.beta;

    num_cout << std::endl << "The value of ground energy is " << H.getGroundEnergy() << std::endl;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    num_cout << "<H> = " << rho.getAverageEnergy() << std::endl;

    std::ofstream out;
    out.open("output/Stat.En.dat");
    out << iomanip_prefs << rho.getAverageEnergy() << std::endl;
    out.close();

    //finishing of creation
    std::cout << std::endl;
    std::cout << "All parts are created!" << std::endl;
    std::cout << std::endl;

    FieldOperatorContainer Operators(S,H,IndexInfo);

    // Single-particle Green's functions
    GFContainer G(S,H,rho,IndexInfo,Operators);
    G.prepare();
    G.compute();
    writeGForHartmut(G,opt.NumberOfMatsubaras);

    /*
    printFramed("Two Particle Green's function calculation");

    std::vector<TwoParticleGFContainer::IndexCombination*> v1;
    v1.push_back(new TwoParticleGFContainer::IndexCombination(0,IndexInfo.getIndexSize()/2,0,IndexInfo.getIndexSize()/2));
    v1.push_back(new TwoParticleGFContainer::IndexCombination(0,0,0,0));
    TwoParticleGFContainer Chi4(S,H,rho,IndexInfo,Operators);
    Chi4.readInitialIndices(v1);
    Chi4.prepare();
    Chi4.compute(wn);
    saveChi("Chi4.dat",Chi4,wn);
    
    Vertex4 Gamma4(IndexInfo,Chi4,G);
    Gamma4.prepareUnAmputated();
    Gamma4.computeUnAmputated();
    Gamma4.prepareAmputated(v1);
    Gamma4.computeAmputated();
    saveGamma("Gamma4.dat",Gamma4,v1,wn,Amputated);*/

    return 0;
};
