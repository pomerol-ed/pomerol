/** \file src/DensityMatrix.cpp
** \brief Density matrix of the grand canonical ensemble.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "DensityMatrix.h"

DensityMatrix::DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta) : 
    S(S), H(H), parts(NULL), beta(beta)
{}

DensityMatrix::~DensityMatrix()
{
    if(parts)
        for(BlockNumber n = 0; n < NumOfBlocks; n++) delete parts[n];
    delete[] parts;
}

void DensityMatrix::prepare(void)
{
    NumOfBlocks = S.NumberOfBlocks();
    RealType GroundEnergy = H.getGroundEnergy();
    parts = new DensityMatrixPart* [NumOfBlocks];
    // There is one-to-one correspondence between parts of the Hamiltonian
    // and parts of the density matrix itself. 
    for(BlockNumber n = 0; n < NumOfBlocks; n++)
        parts[n] = new DensityMatrixPart(H.part(n),beta,GroundEnergy);
}

void DensityMatrix::compute(void)
{
    RealType Z = 0;
    // A total partition function is a sum over partition functions of
    // all non-normalized parts.
    for(BlockNumber n = 0; n < NumOfBlocks; n++) Z += parts[n]->compute();
    // Divide the density matrix by Z.
    for(BlockNumber n = 0; n < NumOfBlocks; n++) parts[n]->normalize(Z);
    cout << "Partition Function = " << Z << endl;
}

RealType DensityMatrix::operator()( QuantumState &state )
{
    BlockNumber block_num = S.getBlockNumber(S.getStateInfo(state));
    int inner_state = S.getInnerState(state);

    return parts[block_num]->weight(inner_state);
}

DensityMatrixPart& DensityMatrix::part(const QuantumNumbers &in)
{
    return *parts[S.getBlockNumber(in)];
}

DensityMatrixPart& DensityMatrix::part(BlockNumber in)
{
    return *parts[in];
}

RealType DensityMatrix::getBeta()
{
    return beta;
}

RealType DensityMatrix::getAverageEnergy()
{
    RealType E = 0;
    for(BlockNumber n = 0; n < NumOfBlocks; n++) E += parts[n]->getAverageEnergy();
    return E;
};
