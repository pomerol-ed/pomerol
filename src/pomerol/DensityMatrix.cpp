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

void DensityMatrix::truncateBlocks(RealType Tolerance, bool verbose)
{
    // init vector<bool> at false
    block_retained.resize(S.NumberOfBlocks());
    for(int i=0; i<block_retained.size(); i++)  block_retained[i]=false;

    // check if each state has weight larger than Tolerance,
    // and set true for the corresponding block
    RealType sum=0;  // sum of weights (to be 1.0)
    int n_states_relevant=0;  // count states having weight larger than Tolerance
    for (QuantumState i=0; i<S.getNumberOfStates(); i++){
        RealType w = getWeight(i);
        if ( w > Tolerance ){
            block_retained[S.getBlockNumber(i)] = true;
            sum += w;
            ++n_states_relevant;
        }
    }

    if(verbose){
        // count blocks and states retained
        int n_blocks_retained=0, n_states_retained=0;
        for(BlockNumber i=0; i<S.NumberOfBlocks(); i++)
            if(isRetained(i)){
                ++n_blocks_retained;
                n_states_retained += S.getBlockSize(i);
            }
        INFO("Number of blocks retained: " << n_blocks_retained);
        INFO("Number of states retained: " << n_states_retained);
        INFO("Number of states having weight larger than Tolerance: " << n_states_relevant);
        INFO("Total weight of the trucated states: " << sum-1.0);
    }
}

bool DensityMatrix::isRetained(BlockNumber in) const
{
    if (int(in) < block_retained.size())  // if function truncateBlocks has been called
        return block_retained[in];
    else
        return true;
}

} // end of namespace Pomerol
