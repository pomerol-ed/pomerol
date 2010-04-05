#include "config.h"
#include "getstates.h"
#include "hpart.h"
#include "CCXpair.h"
#include "output.h"
#include "green.h"
#include "hamiltonian.h"

#include "iniconfig.h"

string input = "system.ini";

getStates S; 
getHpart M(S);
output_handle OUT;
hamiltonian H(S,OUT,input);

//(S,OUT,"system.ini",true,true);
//vector<getHpart> **Hpart;


matrixs operatorsCiCXj(S,H,OUT);


int main()
{	
	//determination of system
	
	IniConfig Ini(input);
	
	int N_bit = Ini["system:N_bit"];
	int N_bit_m = Ini["system:N_bit_m"];
	OUT = output_handle(Ini["output:path"]);

	// parameters of Green Function
	
	int i = Ini["Green Function:i"];
	int j = Ini["Green Function:j"];

	S.inigetStates(N_bit, N_bit_m);

	
	//end of determination	
	
	cout << "=======================" << endl;
	cout << "System is determinated" << endl;
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
	return 0;
}







