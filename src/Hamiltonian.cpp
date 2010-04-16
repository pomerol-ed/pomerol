#include "Hamiltonian.h"
#include "iniconfig.h"

extern IniConfig* pIni;

Hamiltonian::Hamiltonian(BitClassification &F_,StatesClassification &S_,output_handle &OUT_, string& config_path_):Formula(F_),S(S_),OUT(OUT_),config_path(config_path_){}

void Hamiltonian::enter(bool diag, bool dump)
{
	IniConfig Ini(config_path);

	output_handle OUT_EVec(OUT.path()+"/EigenVec");		// create output_folders
	output_handle OUT_EVal(OUT.path()+"/EigenVal");		// create output_folders

	RealType J = (*pIni)["system:J_c"];
	RealType U = (*pIni)["system:U_c"];
	RealType Us = (*pIni)["system:Us_c"];
	RealType mu = (*pIni)["system:mu_c"];
	RealType mus = (*pIni)["system:mus_c"];
	RealType t = (*pIni)["system:t_c"];
	RealType ts = (*pIni)["system:ts_c"];

	Hpart = new HamiltonianPart * [S.NumberOfBlocks()];
	for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
	{
	  Hpart[current_block] = new HamiltonianPart(Formula,S,S.getBlockInfo(current_block));
	  Hpart[current_block]->iniHamiltonianPart( J, U, Us, mu, mus, t, ts, OUT_EVal.path(), OUT_EVec.path());

	  if (diag) Hpart[current_block]->diagonalization();
	  if (dump) Hpart[current_block]->dump();
	  cout << "Hpart" << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is entered";
	  if (diag) cout << ", diagonalized" << flush;
	  if (dump) cout << " and dumped";
	  cout << ". Size = " << S.clstates(S.getBlockInfo(current_block)).size() << endl; 

	}
}

HamiltonianPart& Hamiltonian::part(const QuantumNumbers &in)
{
  return *Hpart[S.getBlockNumber(in)];
}

HamiltonianPart& Hamiltonian::part(BlockNumber in)
{
  return *Hpart[in];
}

RealType Hamiltonian::eigenval(QuantumState &state)
{
 int inner_state = S.inner_state(state);
 return part(S.getStateInfo(state)).reV(inner_state);
}

void Hamiltonian::diagonalize()
{
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
    {
      Hpart[current_block]->diagonalization();
      cout << "Hpart" << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is diagonalized." << endl;
    }
}

void Hamiltonian::dump()
{
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
    {
      Hpart[current_block]->dump();
    }
  cout << "Hamiltonian has been dumped." << endl;
}
