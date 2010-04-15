#include "DensityMatrix.h"

DensityMatrix::DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta) : S(S)
{
    NumOfBlocks = S.NumberOfBlocks();
    parts = new DensityMatrixPart* [NumOfBlocks];
    
    RealType Z = 0;
    for(BlockNumber n = 0; n < NumOfBlocks; n++){
        parts[n] = new DensityMatrixPart(H.block(n),beta);
        Z += parts[n]->getZ();
    }
    for(BlockNumber n = 0; n < NumOfBlocks; n++) parts[n]->normalize(Z);
}

DensityMatrix::~DensityMatrix()
{
    delete[] parts;
}

RealType DensityMatrix::operator()( QuantumState &state )
{
    int inner_state = S.inner_state(state);
    BlockNumber block_num = S.getBlockNumber(S.getStateInfo(state));
    
    return parts[block_num]->weight(inner_state);
}

DensityMatrixPart& DensityMatrix::block(const QuantumNumbers &in)
{
    return *parts[S.getBlockNumber(in)];
}

DensityMatrixPart& DensityMatrix::block(BlockNumber in)
{
    return *parts[in];
}