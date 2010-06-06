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
                AnnihilationOperator& C0, AnnihilationOperator& C1, 
                CreationOperator& CX2, CreationOperator& CX3,
                DensityMatrix& DM,
                output_handle &OUT) : parts(0), S(S), H(H), C0(C0), C1(C1), CX2(CX2), CX3(CX3), DM(DM)
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
    switch(getPermutation3(OperatorPosition).perm[OperatorPosition]){
      case 0: return C0.mapsTo(in);
      case 1: return C1.mapsTo(in);
      case 2: return CX2.mapsTo(in);
      default: return ERROR_BLOCK_NUMBER;
    }
}

void Vertex4::prepare(void)
{
    std::list<BlockMapping> CX3NontrivialBlocks = CX3.getNonTrivialIndices();
  
    for(std::list<BlockMapping>::const_iterator outer_iter = CX3NontrivialBlocks.begin();
        outer_iter != CX3NontrivialBlocks.end(); outer_iter++){
            for(size_t p=0; p<6; ++p){ // Search for non-vanishing world lines
                  BlockNumber blocks[4];
                  blocks[3] = outer_iter->second;
                  blocks[2] = outer_iter->first;
                  blocks[1] = OperatorAtPositionMapsTo(p,2,blocks[2]);
                  blocks[0] = OperatorAtPositionMapsTo(p,1,blocks[1]);
                  if(OperatorAtPositionMapsTo(p,0,blocks[0]) == blocks[4]){
                      //parts.push_back(new Vertex4Part(
                      //    
                      //,p));
                  }
            }
    }  
}

void Vertex4::compute(void)
{
    for(std::list<Vertex4Part*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        (*iter)->compute();
}

ComplexType Vertex4::operator()(ComplexType Frequency0, ComplexType Frequency1, ComplexType Frequency2)
{
    ComplexType Value = 0;
    for(std::list<Vertex4Part*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        Value += (**iter)(Frequency0,Frequency1,Frequency2);
    return Value;
}

string Vertex4::getPath()
{
    return green_path.fullpath();
}