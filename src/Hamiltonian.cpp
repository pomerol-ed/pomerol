#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(BitClassification &F_,StatesClassification &S_,output_handle &OUT_, string& config_path_):Formula(F_),S(S_),OUT(OUT_),config_path(config_path_){}

void Hamiltonian::enter()
{
	output_handle OUT_EVec(OUT.path()+"/EigenVec");		// create output_folders
	output_handle OUT_EVal(OUT.path()+"/EigenVal");		// create output_folders

	Hpart = new HamiltonianPart * [S.NumberOfBlocks()];
	for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
	{
	  Hpart[current_block] = new HamiltonianPart(Formula,S,S.getBlockInfo(current_block),OUT_EVal.path(), OUT_EVec.path());
	  Hpart[current_block]->enter();

	  cout << "Hamiltonian block " << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is entered";
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
 int inner_state = S.getInnerState(state);
 return part(S.getStateInfo(state)).reV(inner_state);
}

void Hamiltonian::diagonalize()
{
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
    {
      Hpart[current_block]->diagonalization();
      cout << "Hpart" << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is diagonalized." << endl;
    }
 computeGroundEnergy();
}

void Hamiltonian::dump()
{
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
    {
      Hpart[current_block]->dump();
    }
  cout << "Hamiltonian has been dumped." << endl;
}

void Hamiltonian::computeGroundEnergy()
{
  RealVectorType LEV(S.NumberOfBlocks());
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
  {
	  LEV(current_block)=Hpart[current_block]->getMinimumEigenvalue();
  }
  GroundEnergy=LEV.minCoeff();
}

RealType Hamiltonian::getGroundEnergy()
{
  return GroundEnergy;
};

void Hamiltonian::reduce(const RealType Cutoff)
{
	cout << "Performing EV cutoff at " << Cutoff << " level" << endl;
	for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
  	{
			Hpart[current_block]->reduce(GroundEnergy+Cutoff);
  	}

};


