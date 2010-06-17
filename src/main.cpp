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
#include "TwoParticleGF.h"

#include "iniconfig.h"


string input = "system.ini";

LatticeAnalysis Lattice;
BitClassification Formula(Lattice);
StatesClassification S(Formula); 
output_handle OUT;
Hamiltonian H(Formula,S,OUT,input);

IniConfig* pIni;

extern long term_counter;

int main()
{	
	cout << "=======================" << endl;
	cout << "Lattice Info" << endl;
	cout << "=======================" << endl;
	Lattice.readin();
	cout << Lattice.printSitesList().str() << flush;
	Formula.prepare();
	cout << "=======================" << endl;
	cout << "System Info" << endl;
	cout << "=======================" << endl;
	Formula.printBitInfoList();
	cout << "=======================" << endl;
	cout << "Hopping Matrix" << endl;
	cout << "=======================" << endl;
	Formula.printHoppingMatrix();
	cout << "=======================" << endl;
	cout << "Terms check" << endl;
	cout << "=======================" << endl;
	Formula.printTerms();
	//determination of system
	
	pIni = new IniConfig(input);
	
	OUT = output_handle((*pIni)["output:path"]);

	// parameters of Green Function
	
	int i = (*pIni)["Green Function:i"];
	int j = (*pIni)["Green Function:j"];

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
	H.diagonalize();
	RealType beta = (*pIni)["Green Function:beta"];
	RealType ProbabilityCutoff = RealType((*pIni)["system:ProbabilityCutoff"]);
	if (ProbabilityCutoff>0 && ProbabilityCutoff < 1) H.reduce(-log(ProbabilityCutoff)/beta);
	H.dump();

	cout << endl << "The value of ground energy is " << H.getGroundEnergy() << endl;

    DensityMatrix rho(S,H,beta);
//	DensityMatrix.reduce();
    	rho.prepare();
    	rho.compute();
	
	//finishing of creation
	cout << endl;
	cout << "All parts are created!" << endl;
	cout << endl;
	
	//begining of creation matrixes C and CX
	
	cout << endl;
	cout << "==========================================" << endl;
	cout << "Beginning of rotation of matrices C and CX" << endl;
	CreationOperator CX(S,H,OUT,i);
	CX.prepare();
	CX.compute();
	CX.dump();
    
	AnnihilationOperator C(S,H,OUT,j);
	C.prepare();
	C.compute();
	C.dump();

	std::list<BlockMapping> ind = CX.getNonTrivialIndices();

	std::list<std::pair<BlockNumber,BlockNumber> >::iterator it;
	for (it=ind.begin();it!=ind.end();it++)
	{ cout << (*it).first << "->" << (*it).second << " = " << CX.getLeftIndex((*it).second) << "->" << CX.getRightIndex((*it).first) << endl;
	}
    // DEBUG
    
	cout << endl;
	cout << "==========================================" << endl;
	cout << "Calculating G_{" << i << j << "}" << endl;
	cout << "==========================================" << endl;
  	GreensFunction G(S,H,C,CX,rho,OUT);
    	G.prepare();
        G.compute();
        
        //std::list<GreensFunctionPart::GreensTerm> terms = G.getTerms();
        //for(std::list<GreensFunctionPart::GreensTerm>::iterator term = terms.begin(); term != terms.end(); term++)
        //    DEBUG(*term)
        
    	G.dumpMatsubara((int)(*pIni)["Green Function:points"]);
    	cout << endl << "All done." << endl;
    
	cout << endl;
	cout << "==========================================" << endl;
	cout << "Two Particle Green's function calculation" << endl;
	cout << "==========================================" << endl;
        TwoParticleGF Chi4(S,H,C,C,CX,CX,rho,OUT);
        Chi4.prepare();
        Chi4.compute(60);

	cout << term_counter << " terms" <<  endl;
	cout << Chi4(0,2,0) << endl;
	cout << Chi4(0,1,5) << endl;

    return 0;
}
