#include "TwoParticleGF.h"

inline Permutation3 getPermutation3(size_t p)
{
    static Permutation3 perms[6] = {
        {{0,1,2},1}, {{0,2,1},-1}, {{1,0,2},-1}, {{1,2,0},1}, {{2,0,1},1}, {{2,1,0},-1}
    };

    return perms[p];
}

/*
inline Permutation4 getPermutation4(size_t p)
{
    static Permutation4 perms[24] = {
        {{0,1,2,3},1}, {{0,1,3,2},-1}, {{0,2,1,3},-1}, {{0,2,3,1},1},
        {{0,3,1,2},1}, {{0,3,2,1},-1}, {{1,0,2,3},-1}, {{1,0,3,2},1},
        {{1,2,0,3},1}, {{1,2,3,0},-1}, {{1,3,0,2},-1}, {{1,3,2,0},1},
        {{2,0,1,3},1}, {{2,0,3,1},-1}, {{2,1,0,3},-1}, {{2,1,3,0},1},
        {{2,3,0,1},1}, {{2,3,1,0},-1}, {{3,0,1,2},-1}, {{3,0,2,1},1},
        {{3,1,0,2},1}, {{3,1,2,0},-1}, {{3,2,0,1},-1}, {{3,2,1,0},1}
    };
  
    return perms[p];
}
*/

extern IniConfig* pIni;

TwoParticleGF::TwoParticleGF(StatesClassification& S, Hamiltonian& H,
                AnnihilationOperator& C1, AnnihilationOperator& C2, 
                CreationOperator& CX3, CreationOperator& CX4,
                DensityMatrix& DM,
                output_handle &OUT) : parts(0), S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM)
{
    green_path = output_handle(OUT.path() + "/Gamma4");
	Storage = new TwoParticleGFPart::MatsubaraContainer(DM.getBeta());
}

TwoParticleGF::~TwoParticleGF()
{
      for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

BlockNumber TwoParticleGF::getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.getLeftIndex(RightIndex);
      case 1: return C2.getLeftIndex(RightIndex);
      case 2: return CX3.getLeftIndex(RightIndex);
      default: return ERROR_BLOCK_NUMBER;
    }
}

BlockNumber TwoParticleGF::getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex)
{
	switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
  		case 0: return C1.getRightIndex(LeftIndex); 
  		case 1: return C2.getRightIndex(LeftIndex);
		case 2: return CX3.getRightIndex(LeftIndex);
  		default: return ERROR_BLOCK_NUMBER;
	}
}

FieldOperatorPart& TwoParticleGF::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.getPartFromLeftIndex(LeftIndex);
      case 1: return C2.getPartFromLeftIndex(LeftIndex);
      case 2: return CX3.getPartFromLeftIndex(LeftIndex);
      default: assert(0);
    }
}

void TwoParticleGF::prepare(void)
{
    std::list<BlockMapping> CX4NontrivialBlocks = CX4.getNonTrivialIndices();
  
    for(std::list<BlockMapping>::const_iterator outer_iter = CX4NontrivialBlocks.begin();
        outer_iter != CX4NontrivialBlocks.end(); outer_iter++){
            for(size_t p=0; p<6; ++p){ // Search for non-vanishing world lines
                  BlockNumber LeftIndices[4];
                  LeftIndices[0] = outer_iter->second;
                  LeftIndices[3] = outer_iter->first;
                  LeftIndices[2] = getLeftIndex(p,2,LeftIndices[3]);
                  LeftIndices[1] = getRightIndex(p,0,LeftIndices[0]);
                  if(getRightIndex(p,1,LeftIndices[1]) == LeftIndices[2] && LeftIndices[1].isCorrect() && LeftIndices[2].isCorrect()){
                      // DEBUG
                      DEBUG("new part: "  << S.getBlockInfo(LeftIndices[0]) << " " 
                                          << S.getBlockInfo(LeftIndices[1]) << " "
                                          << S.getBlockInfo(LeftIndices[2]) << " "
                                          << S.getBlockInfo(LeftIndices[3]) << " "
                      <<"BlockNumbers part: "  << LeftIndices[0] << " " << LeftIndices[1] << " " << LeftIndices[2] << " " << LeftIndices[3]);
					  //DEBUG(OperatorPartAtPosition(p,1,LeftIndices[1]).getRowMajorValue());
                      parts.push_back(new TwoParticleGFPart(
                            OperatorPartAtPosition(p,0,LeftIndices[0]),
                            OperatorPartAtPosition(p,1,LeftIndices[1]),
                            OperatorPartAtPosition(p,2,LeftIndices[2]),
                            (CreationOperatorPart&)CX4.getPartFromLeftIndex(LeftIndices[3]),
                            H.part(LeftIndices[0]), H.part(LeftIndices[1]), H.part(LeftIndices[2]), H.part(LeftIndices[3]),
                            DM.part(LeftIndices[0]), DM.part(LeftIndices[1]), DM.part(LeftIndices[2]), DM.part(LeftIndices[3]),
                      getPermutation3(p)));
                  }
            }
    }  
}

void TwoParticleGF::compute(long NumberOfMatsubaras)
{
	int i=0;
	Storage->prepare(NumberOfMatsubaras);
    for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
    {
	    cout << (int) ((i*100.0)/parts.size()) << "  " <<flush;
        (*iter)->compute(NumberOfMatsubaras);
		*Storage+=(*iter)->getMatsubaraContainer();
		(*iter)->clear();
		++i;
    }
	cout << endl;
}

ComplexType TwoParticleGF::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    ComplexType Value = 0;
//    for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
//        Value += (**iter)(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
    return (*Storage)(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

}

unsigned short TwoParticleGF::getBit(size_t Position) const
{
    switch(Position){
        case 0: return C1.getBit();
        case 1: return C2.getBit();
        case 2: return CX3.getBit();
        case 3: return CX4.getBit();
        default: assert(0);
    }
}

RealType TwoParticleGF::getBeta() const
{
    return DM.getBeta();
}

string TwoParticleGF::getPath()
{
    return green_path.fullpath();
}
