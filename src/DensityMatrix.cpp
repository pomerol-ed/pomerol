#include "DensityMatrix.h"

DensityMatrix::DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta) : S(S)
{
    BlockNumber num = S.NumberOfBlocks();
    parts.reserve(num);
    
    RealType Z = 0;
    for(BlockNumber n = 0; n < num; n++){
        parts.push_back(DensityMatrixPart(H.block(n),beta));
        Z += parts[n].getZ();
    }
    for(BlockNumber n = 0; n < num; n++) parts[n].normalize(Z);
}

RealType DensityMatrix::operator()( QuantumState &state )
{
    int inner_state = S.inner_state(state);
    BlockNumber block_num = S.getBlockNumber(S.getStateInfo(state));
    
    return parts[block_num].weight(inner_state);
}