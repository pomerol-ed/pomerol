#include "config.h"
#include "output.h"
#include "Term.h"
#include "BitClassification.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "FieldOperatorPart.h"
#include "FieldOperator.h"
#include "GreensFunction.h"
#include "Vertex4.h"

#include "iniconfig.h"


string input = "system.ini";

BitClassification Formula;
StatesClassification S(Formula); 
output_handle OUT;
Hamiltonian H(Formula,S,OUT,input);

IniConfig* pIni;

int main()
{	
	Formula.readin();
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
	H.dump();
	cout << endl << "The value of ground energy is " << H.getGroundEnergy() << endl;

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
	RealType beta = (*pIni)["Green Function:beta"];
    	DensityMatrix rho(S,H,beta);
    	rho.prepare();
    	rho.compute();
    
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
	cout << "4th order Vertex calculation" << endl;
	cout << "==========================================" << endl;
        Vertex4 Gamma4(S,H,C,C,CX,CX,rho,OUT);
        Gamma4.prepare();
        Gamma4.compute();
        
	cout << Gamma4(0,1,3) << endl;

    	return 0;
}
