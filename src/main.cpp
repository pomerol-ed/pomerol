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
	exit(0);
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

//	bool dump_ev_ef = Ini["output:dump_ev_ef"];

	H.enter(true,true);
	//finishing of creation
	
	cout << endl;
	cout << "All parts are created!" << endl;
	cout << endl;
	
/*	if (dump_ev_ef)
	{
		cout << "Eigenvectors are placed in " << OUT_EVec.fullpath() << endl;
		cout << "Eigenvalues are placed in " << OUT_EVal.fullpath() << endl;
	}
	cout << endl;*/

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
    // DEBUG
    RealType beta = (*pIni)["Green Function:beta"];
    DensityMatrix rho(S,H,beta);
    rho.compute();
    
    GreensFunction G(S,H,C,CX,rho,OUT);
    G.compute();
    //StatesClassification& _S(S);
    //CreationOperator& _CX(CX);
    //AnnihilationOperator& _C(C);
    //Hamiltonian& _H(H);
    //DensityMatrix& _rho(rho);
    
//	CX.print_to_screen();


//	test.enter();
	//test.print_to_screen();
	/*
	operatorsCiCXj.inimatrixs(i,j);
	
	operatorsCiCXj.putMatrXCX();				//creation matrix UXCXU

	if (i!=j)
	
{	
	operatorsCiCXj.putMatrXC();				//creation matrix UXCU 
	if ((bool) Ini["output:dump_c_cx"])
	{
		operatorsCiCXj.dump();
	
		cout << endl;
		cout << "Notrivial elements of rotated matrices" << endl;
		cout << "are placed in " << operatorsCiCXj.path() << endl;
	}
}	
	cout << endl;

	//finishing of creation

	//begining of creation Green Function G_ij(w)

	cout << endl;
	cout << "====================================" << endl;
	cout << "Start of Greens Function calculation" << endl;
	cout << endl;

	green GF_ij(i,j,S,H,operatorsCiCXj,OUT);
	
	if ((bool) (*pIni)["output:dump_green"])
	{
		GF_ij.dump();
	
  		cout << endl;
		cout << "GF is placed in " << GF_ij.path() << endl;
		cout << endl;
	}

	//finishing of creation
	*/
	cout << endl << "All done." << endl;

	return 0;
}
