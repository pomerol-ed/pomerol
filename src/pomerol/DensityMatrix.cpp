#include "pomerol/DensityMatrix.h"

namespace Pomerol{

DensityMatrix::DensityMatrix(const StatesClassification& S, const Hamiltonian& H, RealType beta) : 
    Thermal(beta), ComputableObject(), S(S), H(H)
{}

DensityMatrix::~DensityMatrix()
{
    for(std::vector<DensityMatrixPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
	delete *iter;
}

void DensityMatrix::prepare(void)
{
    if (Status >= Prepared) return;
    parts = std::vector<DensityMatrixPart*>(S.NumberOfBlocks());
    BlockNumber NumOfBlocks = parts.size();
    RealType GroundEnergy = H.getGroundEnergy();
    // There is one-to-one correspondence between parts of the Hamiltonian
    // and parts of the density matrix itself. 
    for(BlockNumber n = 0; n < NumOfBlocks; n++)
        parts[n] = new DensityMatrixPart(S, H.getPart(n),beta,GroundEnergy);
    Status = Prepared;
}

void DensityMatrix::compute(void)
{
    if (Status >= Computed) return;
    RealType Z = 0;
    // A total partition function is a sum over partition functions of
    // all non-normalized parts.
    for(std::vector<DensityMatrixPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        Z += (*iter)->computeUnnormalized();
 
    // Divide the density matrix by Z.
    for(std::vector<DensityMatrixPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        (*iter)->normalize(Z);
    Status = Computed;
}

RealType DensityMatrix::getWeight(QuantumState state) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    BlockNumber BlockNumber = S.getBlockNumber(state);
    InnerQuantumState InnerState = S.getInnerState(state);

    return parts[BlockNumber]->getWeight(InnerState);
}
 
const DensityMatrixPart& DensityMatrix::getPart(const QuantumNumbers &in) const
{
    return *parts[S.getBlockNumber(in)];
}
 
const DensityMatrixPart& DensityMatrix::getPart(BlockNumber in) const
{
    return *parts[in];
}

RealType DensityMatrix::getAverageEnergy() const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType E = 0;
    for(std::vector<DensityMatrixPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
    E += (*iter)->getAverageEnergy();
    return E;
};

RealType DensityMatrix::getAverageOccupancy() const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType n = 0;
    for(std::vector<DensityMatrixPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
    n += (*iter)->getAverageOccupancy();
    return n;
};

RealType DensityMatrix::getAverageOccupancy(ParticleIndex i) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType n = 0;
    for(std::vector<DensityMatrixPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
    n += (*iter)->getAverageOccupancy(i);
    return n;
};

RealType DensityMatrix::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType NN = 0;
    for(std::vector<DensityMatrixPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
    NN += (*iter)->getAverageDoubleOccupancy(i,j);
    return NN;
};
 
} // end of namespace Pomerol
