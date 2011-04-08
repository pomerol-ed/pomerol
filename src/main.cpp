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

#include<fstream>
#include "iniconfig.h"

string input = "system.ini";

LatticeAnalysis Lattice;
BitClassification IndexInfo(Lattice);
StatesClassification S(IndexInfo); 
output_handle OUT;
Hamiltonian H(IndexInfo,S,OUT,input);

ostream &OUTPUT_STREAM=std::cout;

IniConfig* pIni;

extern long term_counter;


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
//
        TwoParticleGFContainer Chi4(S,H,rho,IndexInfo,Operators);
        Chi4.readNonTrivialIndices(v1);
        Chi4.prepare();
        Chi4.compute(30);

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
               };

    cout << term_counter << " terms" <<  endl;
    exit(0);
    int i = 1; //(*pIni)["Green Function:i"];
    int j = 1; //(*pIni)["Green Function:j"];
    cout << endl;
    cout << "==========================================" << endl;
    cout << "Beginning of rotation of matrices C and CX" << endl;
    CreationOperator CX(S,H,i);
    CX.prepare();
    CX.compute();
    //CX.dump();
    
    AnnihilationOperator C(S,H,j);
    C.prepare();
    C.compute();
    //C.dump();

    // DEBUG
    /*std::list<BlockMapping> ind = CX.getNonTrivialIndices();

    std::list<std::pair<BlockNumber,BlockNumber> >::iterator it;
    for (it=ind.begin();it!=ind.end();it++)
    { cout << (*it).first << "," << (*it).second << " = " << CX.getLeftIndex((*it).second) << "," << CX.getRightIndex((*it).first) << endl;
    }
*/

/*    cout << endl;
    cout << "==========================================" << endl;
    cout << "Calculating G_{" << i << j << "}" << endl;
    cout << "==========================================" << endl;
      GreensFunction G(S,H,C,CX,rho);
        G.prepare();
        cout << G.vanishes() << endl;
        G.compute();
        
        //std::list<GreensFunctionPart::GreensTerm> terms = G.getTerms();
        //for(std::list<GreensFunctionPart::GreensTerm>::iterator term = terms.begin(); term != terms.end(); term++)
        //    DEBUG(*term)
        
    G.dumpMatsubara((int)(*pIni)["Green Function:points"]);
    cout << endl << "All done." << endl;

    cout << i << j << ": G.vanishes() = " << G.vanishes() << endl;

  */
return 0;
};
/*
    //begining of creation matrixes C and CX
    
    // parameters of Green Function
    
    
    if ((*pIni)["System:calculate_2PGF"]){
        cout << endl;
        cout << "==========================================" << endl;
        cout << "Two Particle Green's function calculation" << endl;
        cout << "==========================================" << endl;

        AnnihilationOperator C1(S,H,OUT,i);
        C1.prepare();
        C1.compute();
        AnnihilationOperator C2(S,H,OUT,j);
        C2.prepare();
        C2.compute();
        CreationOperator CX3(S,H,OUT,i);
        CX3.prepare();
        CX3.compute();
        CreationOperator CX4(S,H,OUT,j);
        CX4.prepare();
        CX4.compute();

        TwoParticleGF CurrentChi4(S,H,C1,C2,CX3,CX4,rho);
        CurrentChi4.prepare();
        CurrentChi4.compute(30);
        //Chi4.dump();

        cout << CurrentChi4(3,2,0) << endl;
        cout << CurrentChi4(2,5,2) << endl;
        cout << CurrentChi4(5,2,2) << endl;
        cout << CurrentChi4(1,7,1) << endl;
        cout << CurrentChi4(2,-2,4) << endl;
        cout << CurrentChi4(29,-29,29) << endl;
        
        

        //Vertex4 Gamma4(Chi4,G,G,G,G);
        //cout << Gamma4(0,2,0) << endl;
        }
    return 0;
}
*/
