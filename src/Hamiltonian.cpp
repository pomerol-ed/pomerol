#include "Hamiltonian.h"
#include "iniconfig.h"

Hamiltonian::Hamiltonian(BitClassification &F_,getStates &S_,output_handle &OUT_, string& config_path_):Formula(F_),S(S_),OUT(OUT_),config_path(config_path_){}

void Hamiltonian::enter(bool diag, bool dump)
{
	IniConfig Ini(config_path);

	output_handle OUT_EVec(OUT.path()+"/EigenVec");		// create output_folders
	output_handle OUT_EVal(OUT.path()+"/EigenVal");		// create output_folders

	RealType J = Ini["system:J_c"];
	RealType U = Ini["system:U_c"];
	RealType Us = Ini["system:Us_c"];
	RealType mu = Ini["system:mu_c"];
	RealType mus = Ini["system:mus_c"];
	RealType t = Ini["system:t_c"];
	RealType ts = Ini["system:ts_c"];

	Hpart = new getHpart * [S.NumberOfBlocks()];
	for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
	{
	  Hpart[current_block] = new getHpart(Formula,S,S.getBlockInfo(current_block));
	  Hpart[current_block]->inigetHpart( J, U, Us, mu, mus, t, ts, OUT_EVal.path(), OUT_EVec.path());

	  if (diag) Hpart[current_block]->diagonalization();
	  if (dump) Hpart[current_block]->dump();
	  cout << "Hpart" << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is entered";
	  if (diag) cout << ", diagonalized" << flush;
	  if (dump) cout << " and dumped";
	  cout << ". Size = " << S.clstates(S.getBlockInfo(current_block)).size() << endl; 

	}
}

getHpart& Hamiltonian::block(const QuantumNumbers &in)
{
  return *Hpart[S.getBlockNumber(in)];
}

getHpart& Hamiltonian::block(BlockNumber in)
{
  return *Hpart[in];
}

RealType Hamiltonian::eigenval(QuantumState &state)
{
 int inner_state = S.inner_state(state);
 return (*this).block(S.getStateInfo(state)).reV(inner_state);
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
