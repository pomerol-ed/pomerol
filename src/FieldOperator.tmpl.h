template<class PartType> 
OperatorContainer<PartType>::OperatorContainer(StatesClassification &System, Hamiltonian &H, output_handle &OUT, int bit) : 
    System(System), H(H), OUT(OUT), bit(bit)
{
    size=0;
}

template<class PartType>
void OperatorContainer<PartType>::prepare()
{
  for (BlockNumber RightIndex=0;RightIndex<System.NumberOfBlocks();RightIndex++)
    {
      BlockNumber LeftIndex = where(RightIndex);
      if (LeftIndex.isCorrect()) 
      {
     	 PartType *Part = new PartType(bit,System,H.part(RightIndex),H.part(LeftIndex),OUT);
	 Data.push_back(Part);
         cout << "Entering " << operatorName << " Operator part " << System.getBlockInfo(RightIndex) << "->" << System.getBlockInfo(LeftIndex) << endl; 
         mapPartsFromRight[RightIndex]=size;
         mapPartsFromLeft[LeftIndex]=size;
	 RightLeftIndices.push_back(std::pair<BlockNumber,BlockNumber>(RightIndex,LeftIndex));
	 LeftRightIndices.push_back(std::pair<BlockNumber,BlockNumber>(LeftIndex,RightIndex));
         size++;
      }    
    }
}

template<class PartType>
std::list<std::pair<BlockNumber,BlockNumber> >& OperatorContainer<PartType>::getNonTrivialIndices(bool OrderFromRight)
{
	return OrderFromRight?RightLeftIndices:LeftRightIndices;
};

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromRightIndex(BlockNumber in)
{
  return *Data[mapPartsFromRight[in]];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromRightIndex(QuantumNumbers in)
{
  return *Data[mapPartsFromRight[System.getBlockNumber(in)]];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromLeftIndex(BlockNumber in)
{
  return *Data[mapPartsFromLeft[in]];
}

template<class PartType>
PartType& OperatorContainer<PartType>::getPartFromLeftIndex(QuantumNumbers in)
{
  return *Data[mapPartsFromLeft[System.getBlockNumber(in)]];
}

template<class PartType>
void OperatorContainer<PartType>::print_to_screen()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
        Data[b_in]->print_to_screen();
  };
}

template<class PartType>
void OperatorContainer<PartType>::compute()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
	cout << (int) ((1.0*b_in/(*this).size) * 100 ) << "  " << flush;
        Data[b_in]->compute();
  };
  cout << endl;
}

template<class PartType>
void OperatorContainer<PartType>::dump()
{
  for (unsigned int b_in=0;b_in<(*this).size;b_in++)
  {
        Data[b_in]->dump();
  };
}
