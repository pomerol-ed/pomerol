#include "Vertex4.h"

inline Permutation3 getPermutation3(size_t p)
{
    static Permutation3 perms[6] = {
        {{0,1,2},1}, {{0,2,1},-1}, {{1,0,2},-1}, {{1,2,0},1}, {{2,0,1},1}, {{2,1,0},-1}
    };
    
    return perms[p];
}

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

extern IniConfig* pIni;

Vertex4::Vertex4(StatesClassification& S, Hamiltonian& H,
                AnnihilationOperator& C1, AnnihilationOperator& C2, 
                CreationOperator& CX3, CreationOperator& CX4,
                DensityMatrix& DM,
                output_handle &OUT) : parts(0), S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM)
{
    green_path = output_handle(OUT.path() + "/Gamma4");
}

Vertex4::~Vertex4()
{
      for(std::list<Vertex4Part*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

BlockNumber Vertex4::OperatorAtPositionMapsTo(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.mapsTo(in); // bad - better use getLeftRightIndices method
      case 1: return C2.mapsTo(in);
      case 2: return CX3.mapsTo(in);
      default: return ERROR_BLOCK_NUMBER;
    }
}

FieldOperatorPart& Vertex4::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber in)
{
    switch(getPermutation3(PermutationNumber).perm[OperatorPosition]){
      case 0: return C1.getPartFromRightIndex(in);
      case 1: return C2.getPartFromRightIndex(in);
      case 2: return CX3.getPartFromRightIndex(in);
      default: assert(0);
    }
}

void Vertex4::prepare(void)
{
    std::list<BlockMapping> CX4NontrivialBlocks = CX4.getNonTrivialIndices();
  
    for(std::list<BlockMapping>::const_iterator outer_iter = CX4NontrivialBlocks.begin();
        outer_iter != CX4NontrivialBlocks.end(); outer_iter++){
            for(size_t p=0; p<6; ++p){ // Search for non-vanishing world lines
                  BlockNumber blocks[4];
                  blocks[3] = outer_iter->second;
                  blocks[2] = outer_iter->first;
                  blocks[1] = OperatorAtPositionMapsTo(p,2,blocks[2]);
                  blocks[0] = OperatorAtPositionMapsTo(p,1,blocks[1]);
                  if(OperatorAtPositionMapsTo(p,0,blocks[0]) == blocks[3]){
                      // DEBUG
                      DEBUG("new part: " << S.getBlockInfo(blocks[0]) << " " << S.getBlockInfo(blocks[1]) << " "<< S.getBlockInfo(blocks[2]) << " "<< S.getBlockInfo(blocks[3]) << " ")
                      parts.push_back(new Vertex4Part(
                            OperatorPartAtPosition(p,0,blocks[0]),
                            OperatorPartAtPosition(p,1,blocks[1]),
                            OperatorPartAtPosition(p,2,blocks[2]),
                            (CreationOperatorPart&)CX4.getPartFromRightIndex(blocks[3]),
                            H.part(blocks[0]), H.part(blocks[1]), H.part(blocks[2]), H.part(blocks[3]),
                            DM.part(blocks[0]), DM.part(blocks[1]), DM.part(blocks[2]), DM.part(blocks[3]),
                      getPermutation3(p)));
                  }
            }
    }  
}

void Vertex4::compute(void)
{
    int i=0;
    for(std::list<Vertex4Part*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
    {
	    DEBUG( i << " / " << parts.size() );
        (*iter)->compute();
	++i;
    }
}

ComplexType Vertex4::operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3)
{
    ComplexType Value = 0;
    for(std::list<Vertex4Part*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        Value += (**iter)(Frequency1,Frequency2,Frequency3);
    return Value;
}

string Vertex4::getPath()
{
    return green_path.fullpath();
}
