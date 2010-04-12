#include "FieldOperator.h"

void CreationOperator::prepare()
{
  Data = new FieldOperatorPart * [System.NumberOfBlocks()];
  for (BlockNumber b=0;b<System.NumberOfBlocks();b++)
    {
      if ((*this).where(b).isCorrect()) 
      {
	 Data[b]=new CreationOperatorPart(bit,System,H.block(b),H.block((*this).where(b)),OUT);
         cout << "Entering Creation Operator part " << System.getBlockInfo(b) << "->" << System.getBlockInfo((*this).where(b)) << endl; 
	 mapNontrivialBlocks[size]=b;
         size++;
      }    
    }
}

FieldOperatorPart& OperatorContainer::part(BlockNumber in)
{
  return *Data[in];
}

FieldOperatorPart& OperatorContainer::part(QuantumNumbers in)
{
  return *Data[System.getBlockNumber(in)];
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

void OperatorContainer::print_to_screen()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
	    Data[mapNontrivialBlocks[b_in]]->print_to_screen();
  };
}

void OperatorContainer::compute()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
	    Data[mapNontrivialBlocks[b_in]]->compute();
  };
}

void OperatorContainer::dump()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
	    Data[mapNontrivialBlocks[b_in]]->dump();
  };
}




