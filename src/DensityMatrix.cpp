/** \file src/DensityMatrix.cpp
** \brief Density matrix of the grand canonical ensemble.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "DensityMatrix.h"

DensityMatrix::DensityMatrix(StatesClassification& S, Hamiltonian& H, RealType beta) : 
    ComputableObject(), S(S), H(H), parts(S.NumberOfBlocks()), beta(beta)
{}

DensityMatrix::~DensityMatrix()
{
    BlockNumber NumOfBlocks = parts.size();
    for(BlockNumber n = 0; n < NumOfBlocks; n++) delete parts[n];
}

void DensityMatrix::prepare(void)
{
    BlockNumber NumOfBlocks = parts.size();
    RealType GroundEnergy = H.getGroundEnergy();
    // There is one-to-one correspondence between parts of the Hamiltonian
    // and parts of the density matrix itself. 
    for(BlockNumber n = 0; n < NumOfBlocks; n++)
        parts[n] = new DensityMatrixPart(S, H.part(n),beta,GroundEnergy);
}

void DensityMatrix::compute(void)
{
    BlockNumber NumOfBlocks = parts.size();
    RealType Z = 0;
    // A total partition function is a sum over partition functions of
    // all non-normalized parts.
    for(BlockNumber n = 0; n < NumOfBlocks; n++) Z += parts[n]->compute();
    // Divide the density matrix by Z.
    for(BlockNumber n = 0; n < NumOfBlocks; n++) parts[n]->normalize(Z);
    std::cout << "Partition Function = " << Z << std::endl;
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

RealType DensityMatrix::getBeta() const
{
    return beta;
}

RealType DensityMatrix::getAverageEnergy()
{
    BlockNumber NumOfBlocks = parts.size();
    RealType E = 0;
    for(BlockNumber n = 0; n < NumOfBlocks; n++) E += parts[n]->getAverageEnergy();
    return E;
};

RealType DensityMatrix::getAverageDoubleOccupancy(ParticleIndex i, ParticleIndex j)
{
    BlockNumber NumOfBlocks = parts.size();
    RealType NN = 0;
    for(BlockNumber n = 0; n < NumOfBlocks; n++) NN += parts[n]->getAverageDoubleOccupancy(i,j);
    return NN;
};


#ifdef pomerolHDF5
void DensityMatrix::save(H5::CommonFG* FG) const
{
    H5::Group RootGroup(FG->createGroup("DensityMatrix"));

    // Dump inverse temperature
    HDF5Storage::saveReal(RootGroup,"beta",beta);

    // Dump parts
    BlockNumber NumOfBlocks = parts.size();
    H5::Group PartsGroup = RootGroup.createGroup("parts");
    for(BlockNumber n = 0; n < NumOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.createGroup(nStr.str().c_str());
	parts[n]->save(&PartGroup);
    }
}

void DensityMatrix::load(const H5::CommonFG* FG)
{
    H5::Group RootGroup(FG->openGroup("DensityMatrix"));  
    beta = HDF5Storage::loadReal(RootGroup,"beta");

    if(Status!=Prepared) prepare();

    H5::Group PartsGroup = RootGroup.openGroup("parts");
    BlockNumber NumOfBlocks = parts.size();
    if(NumOfBlocks != PartsGroup.getNumObjs()){
	// TODO: Something is very wrong here!
	// Throw an exception?
    }
    for(BlockNumber n = 0; n < NumOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.openGroup(nStr.str().c_str());
	parts[n]->load(&PartGroup);
    }
}
#endif // endif :: #ifdef pomerolHDF5
