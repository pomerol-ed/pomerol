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

BlockNumber TwoParticleGF::getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.getLeftIndex(in); // bad - better use getLeftRightIndices method
      case 1: return C2.getLeftIndex(in);
      case 2: return CX3.getLeftIndex(in);
      default: return ERROR_BLOCK_NUMBER;
    }
}

BlockNumber TwoParticleGF::getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in)
{
	switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
  		case 0: return C1.getRightIndex(in); // bad - better use getLeftRightIndices method
  		case 1: return C2.getRightIndex(in);
		case 2: return CX3.getRightIndex(in);
  		default: return ERROR_BLOCK_NUMBER;
	}
}

FieldOperatorPart& TwoParticleGF::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.getPartFromLeftIndex(in);
      case 1: return C2.getPartFromLeftIndex(in);
      case 2: return CX3.getPartFromLeftIndex(in);
      default: assert(0);
    }
}

void TwoParticleGF::prepare(void)
{
    std::list<BlockMapping> CX4NontrivialBlocks = CX4.getNonTrivialIndices();
  
    for(std::list<BlockMapping>::const_iterator outer_iter = CX4NontrivialBlocks.begin();
        outer_iter != CX4NontrivialBlocks.end(); outer_iter++){
            for(size_t p=0; p<6; ++p){ // Search for non-vanishing world lines
                  BlockNumber blocks[4];
                  blocks[0] = outer_iter->second;
                  blocks[3] = outer_iter->first;
                  blocks[2] = getLeftIndex(p,2,blocks[3]);
                  blocks[1] = getRightIndex(p,1,blocks[0]);
                  if(getLeftIndex(p,1,blocks[1]) == blocks[2] && blocks[3].isCorrect() && blocks[2].isCorrect()){
                      // DEBUG
                      DEBUG("new part: "  << S.getBlockInfo(blocks[0]) << " " 
                                          << S.getBlockInfo(blocks[1]) << " "
                                          << S.getBlockInfo(blocks[2]) << " "
                                          << S.getBlockInfo(blocks[3]) << " "
                      <<"BlockNumbers part: "  << blocks[0] << " " << blocks[1] << " " << blocks[2] << " " << blocks[3]);
                      parts.push_back(new TwoParticleGFPart(
                            OperatorPartAtPosition(p,0,blocks[0]),
                            OperatorPartAtPosition(p,1,blocks[1]),
                            OperatorPartAtPosition(p,2,blocks[2]),
                            (CreationOperatorPart&)CX4.getPartFromLeftIndex(blocks[3]),
                            H.part(blocks[0]), H.part(blocks[1]), H.part(blocks[2]), H.part(blocks[3]),
                            DM.part(blocks[0]), DM.part(blocks[1]), DM.part(blocks[2]), DM.part(blocks[3]),
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
