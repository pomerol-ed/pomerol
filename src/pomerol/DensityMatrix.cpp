#include "pomerol/DensityMatrix.h"

namespace Pomerol{

template<bool Complex>
DensityMatrix<Complex>::DensityMatrix(const StatesClassification<Complex>& S,
                                      const Hamiltonian<Complex>& H,
                                      RealType beta) :
    Thermal(beta), ComputableObject(), S(S), H(H)
{}

template<bool Complex>
DensityMatrix<Complex>::~DensityMatrix()
{
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
	delete *iter;
}

template<bool Complex>
void DensityMatrix<Complex>::prepare(void)
{
    if (Status >= Prepared) return;
    parts = std::vector<PartT*>(S.NumberOfBlocks());
    BlockNumber NumOfBlocks = parts.size();
    RealType GroundEnergy = H.getGroundEnergy();
    // There is one-to-one correspondence between parts of the Hamiltonian
    // and parts of the density matrix itself.
    for(BlockNumber n = 0; n < NumOfBlocks; n++)
        parts[n] = new PartT(S, H.getPart(n),beta,GroundEnergy);
    Status = Prepared;
}

template<bool Complex>
void DensityMatrix<Complex>::compute(void)
{
    if (Status >= Computed) return;
    RealType Z = 0;
    // A total partition function is a sum over partition functions of
    // all non-normalized parts.
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
        Z += (*iter)->computeUnnormalized();

    // Divide the density matrix by Z.
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
        (*iter)->normalize(Z);
    Status = Computed;
}

template<bool Complex>
RealType DensityMatrix<Complex>::getWeight(QuantumState state) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    BlockNumber BlockNumber = S.getBlockNumber(state);
    InnerQuantumState InnerState = S.getInnerState(state);

    return parts[BlockNumber]->getWeight(InnerState);
}

template<bool Complex>
auto DensityMatrix<Complex>::getPart(const QuantumNumbers<Complex> &in) const -> const PartT&
{
    return *parts[S.getBlockNumber(in)];
}

template<bool Complex>
auto DensityMatrix<Complex>::getPart(BlockNumber in) const -> const PartT&
{
    return *parts[in];
}

template<bool Complex>
RealType DensityMatrix<Complex>::getAverageEnergy() const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType E = 0;
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
    E += (*iter)->getAverageEnergy();
    return E;
};

template<bool Complex>
RealType DensityMatrix<Complex>::getAverageOccupancy() const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType n = 0;
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
    n += (*iter)->getAverageOccupancy();
    return n;
};

template<bool Complex>
RealType DensityMatrix<Complex>::getAverageOccupancy(ParticleIndex i) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType n = 0;
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
    n += (*iter)->getAverageOccupancy(i);
    return n;
};

template<bool Complex>
RealType DensityMatrix<Complex>::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j) const
{
    if ( Status < Computed ) { ERROR("DensityMatrix is not computed yet."); throw (exStatusMismatch()); };
    RealType NN = 0;
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
    NN += (*iter)->getAverageDoubleOccupancy(i,j);
    return NN;
};

template<bool Complex>
void DensityMatrix<Complex>::truncateBlocks(RealType Tolerance, bool verbose)
{
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
        (*iter)->truncate(Tolerance);

    if(verbose){
        // count retained blocks and states included in those blocks
        int n_blocks_retained=0, n_states_retained=0;
        for(BlockNumber i=0; i<S.NumberOfBlocks(); i++)
            if(isRetained(i)){
                ++n_blocks_retained;
                n_states_retained += S.getBlockSize(i);
            }
        INFO("Number of blocks retained: " << n_blocks_retained);
        INFO("Number of states retained: " << n_states_retained);
    }
}

template<bool Complex>
bool DensityMatrix<Complex>::isRetained(BlockNumber in) const
{
    return parts[in]->isRetained();
}

template class DensityMatrix<false>;
template class DensityMatrix<true>;

} // end of namespace Pomerol
