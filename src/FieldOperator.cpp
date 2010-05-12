#include "FieldOperator.h"

CreationOperator::CreationOperator(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit) : 
    OperatorContainer<CreationOperatorPart>(System,H,OUT,bit)
{
    operatorName = "Creation";
}

AnnihilationOperator::AnnihilationOperator(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit) : 
    OperatorContainer<AnnihilationOperatorPart>(System,H,OUT,bit)
{
    operatorName = "Annihilation";
}

QuantumNumbers CreationOperator::where(QuantumNumbers in) // Require explicit knowledge of QuantumNumbers structure - Not very good
{
  int lz, spin;
  System.getSiteInfo(bit,lz,spin);
  QuantumNumbers q_out;
  if (spin == 1) 
    q_out = QuantumNumbers(in.Lz + lz,in.N_up+1,in.N_down);
  else 
    q_out = QuantumNumbers(in.Lz + lz,in.N_up,in.N_down+1);
  return q_out;
}

BlockNumber CreationOperator::where(BlockNumber in)
{
  QuantumNumbers q_in = System.getBlockInfo(in);	
  QuantumNumbers q_out = (*this).where(q_in);
  BlockNumber out = System.getBlockNumber(q_out);
  return out;
}

QuantumNumbers AnnihilationOperator::where(QuantumNumbers in) // Require explicit knowledge of QuantumNumbers structure - Not very good
{
  int lz, spin;
  System.getSiteInfo(bit,lz,spin);
  QuantumNumbers q_out;
  if (spin == 1) 
    q_out = QuantumNumbers(in.Lz - lz,in.N_up-1,in.N_down);
  else 
    q_out = QuantumNumbers(in.Lz - lz,in.N_up,in.N_down-1);
  return q_out;
}

BlockNumber AnnihilationOperator::where(BlockNumber in)
{
  QuantumNumbers q_in = System.getBlockInfo(in);	
  QuantumNumbers q_out = (*this).where(q_in);
  BlockNumber out = System.getBlockNumber(q_out);
  return out;
}
