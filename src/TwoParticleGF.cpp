/** \file src/TwoParticleGF.cpp
** \brief Two-particle Green's function in the Matsubara representation.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "TwoParticleGF.h"
extern output_handle OUT;
extern std::ostream& OUTPUT_STREAM;

static const Permutation3 permutations3[6] = {
    {{0,1,2},1},
    {{0,2,1},-1},
    {{1,0,2},-1},
    {{1,2,0},1},
    {{2,0,1},1},
    {{2,1,0},-1}
};

TwoParticleGF::TwoParticleGF(StatesClassification& S, Hamiltonian& H,
                AnnihilationOperator& C1, AnnihilationOperator& C2, 
                CreationOperator& CX3, CreationOperator& CX4,
                DensityMatrix& DM) : S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM), parts(0)
{
    green_path = output_handle(OUT.path() + "/Gamma4");
    Storage = new TwoParticleGFPart::MatsubaraContainer(DM.getBeta());
    vanish = true;
    Status = Constructed;
}

TwoParticleGF::~TwoParticleGF()
{
      for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

BlockNumber TwoParticleGF::getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex)
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getLeftIndex(RightIndex);
        case 1: return C2.getLeftIndex(RightIndex);
        case 2: return CX3.getLeftIndex(RightIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

BlockNumber TwoParticleGF::getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex)
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getRightIndex(LeftIndex);
        case 1: return C2.getRightIndex(LeftIndex);
        case 2: return CX3.getRightIndex(LeftIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

FieldOperatorPart& TwoParticleGF::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex)
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getPartFromLeftIndex(LeftIndex);
        case 1: return C2.getPartFromLeftIndex(LeftIndex);
        case 2: return CX3.getPartFromLeftIndex(LeftIndex);
        default: assert(0);
    }
}

void TwoParticleGF::prepare(void)
{
if (Status < Prepared)
{
    // Find out non-trivial blocks of CX4.
    std::list<BlockMapping> CX4NontrivialBlocks = CX4.getNonTrivialIndices();
    for(std::list<BlockMapping>::const_iterator outer_iter = CX4NontrivialBlocks.begin();
        outer_iter != CX4NontrivialBlocks.end(); outer_iter++){ // Iterate over the outermost index.
            for(size_t p=0; p<6; ++p){ // Choose a permutation
                  BlockNumber LeftIndices[4];
                  LeftIndices[0] = outer_iter->second;
                  LeftIndices[3] = outer_iter->first;
                  LeftIndices[2] = getLeftIndex(p,2,LeftIndices[3]);
                  LeftIndices[1] = getRightIndex(p,0,LeftIndices[0]);
                  // < LeftIndices[0] | O_1 | LeftIndices[1] >
                  // < LeftIndices[1] | O_2 | getRightIndex(p,1,LeftIndices[1]) >
                  // < LeftIndices[2]| O_3 | LeftIndices[3] >
                  // < LeftIndices[3] | CX4 | LeftIndices[0] >
                  // Select a relevant 'world stripe' (sequence of blocks).
                  if(getRightIndex(p,1,LeftIndices[1]) == LeftIndices[2] && LeftIndices[1].isCorrect() && LeftIndices[2].isCorrect()){
                      // DEBUG
                      DEBUG("new part: "  << S.getBlockInfo(LeftIndices[0]) << " " 
                                          << S.getBlockInfo(LeftIndices[1]) << " "
                                          << S.getBlockInfo(LeftIndices[2]) << " "
                                          << S.getBlockInfo(LeftIndices[3]) << " "
                      <<"BlockNumbers part: "  << LeftIndices[0] << " " << LeftIndices[1] << " " << LeftIndices[2] << " " << LeftIndices[3]);
                      parts.push_back(new TwoParticleGFPart(
                            OperatorPartAtPosition(p,0,LeftIndices[0]),
                            OperatorPartAtPosition(p,1,LeftIndices[1]),
                            OperatorPartAtPosition(p,2,LeftIndices[2]),
                            (CreationOperatorPart&)CX4.getPartFromLeftIndex(LeftIndices[3]),
                            H.part(LeftIndices[0]), H.part(LeftIndices[1]), H.part(LeftIndices[2]), H.part(LeftIndices[3]),
                            DM.part(LeftIndices[0]), DM.part(LeftIndices[1]), DM.part(LeftIndices[2]), DM.part(LeftIndices[3]),
                      permutations3[p]));
                      }
            }
    }  
    if ( parts.size() > 0 ) vanish = false;
    Status = Prepared;
};
}

bool TwoParticleGF::vanishes()
{
return vanish;
}

void TwoParticleGF::compute(long NumberOfMatsubaras)
{
if (Status < Computed){
    Storage->prepare(NumberOfMatsubaras);
    for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
    {
        // TODO: More elegant output.
        std::cout << static_cast<int>((distance(parts.begin(),iter)*100.0)/parts.size()) << "  " << std::flush;
        (*iter)->compute(NumberOfMatsubaras);
        *Storage+=(*iter)->getMatsubaraContainer();
        (*iter)->clear();
    }
    std::cout << std::endl;
    Status = Computed;
    }
}

size_t TwoParticleGF::getNumResonantTerms() const
{
    size_t num = 0;
    for(std::list<TwoParticleGFPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
        num += (*iter)->getNumResonantTerms();
    return num;
}

size_t TwoParticleGF::getNumNonResonantTerms() const
{
    size_t num = 0;
    for(std::list<TwoParticleGFPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
        num += (*iter)->getNumNonResonantTerms();
    return num;
}

ComplexType TwoParticleGF::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    ComplexType Value = 0;
    for(std::list<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++){
        Value+=(**iter)(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3 );
        };
    return (*Storage)(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3)+Value;

}

unsigned short TwoParticleGF::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C1.getIndex();
        case 1: return C2.getIndex();
        case 2: return CX3.getIndex();
        case 3: return CX4.getIndex();
        default: assert(0);
    }
}

RealType TwoParticleGF::getBeta() const
{
    return DM.getBeta();
}

std::string TwoParticleGF::getPath()
{
    return green_path.fullpath();
}
