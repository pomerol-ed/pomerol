#include "Hamiltonian.h"

#warning TODO: Hamiltonian still needs prepare() and compute() methods.

Hamiltonian::Hamiltonian(IndexClassification &F_,StatesClassification &S_,output_handle &OUT_, std::string& config_path_):
ComputableObject(),Formula(F_),S(S_),OUT(OUT_){}

Hamiltonian::~Hamiltonian()
{
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber current_block=0;current_block<NumberOfBlocks;current_block++)
	delete parts[current_block];
}

void Hamiltonian::enter()
{
    output_handle OUT_EVec(OUT.path()+"/EigenVec");		// create output_folders
    output_handle OUT_EVal(OUT.path()+"/EigenVal");		// create output_folders

    BlockNumber NumberOfBlocks = S.NumberOfBlocks();
    parts.resize(NumberOfBlocks);
    for (BlockNumber current_block=0;current_block<NumberOfBlocks;current_block++)
    {
	parts[current_block] = new HamiltonianPart(Formula,S,S.getBlockInfo(current_block),OUT_EVal.path(), OUT_EVec.path());
	parts[current_block]->enter();

	std::cout << "Hamiltonian block " << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is entered";
	std::cout << ". Size = " << S.clstates(S.getBlockInfo(current_block)).size() << std::endl;
    }
    Status = Prepared;
}

HamiltonianPart& Hamiltonian::part(const QuantumNumbers &in)
{
  return *parts[S.getBlockNumber(in)];
}

HamiltonianPart& Hamiltonian::part(BlockNumber in)
{
  return *parts[in];
}

RealType Hamiltonian::eigenval(QuantumState &state)
{
    int inner_state = S.getInnerState(state);
    return part(S.getStateInfo(state)).reV(inner_state);
}

void Hamiltonian::diagonalize()
{
  BlockNumber NumberOfBlocks = parts.size();
  for (BlockNumber current_block=0;current_block<NumberOfBlocks;current_block++)
    {
      parts[current_block]->diagonalization();
      std::cout << "Hpart" << S.getBlockInfo(current_block) << " ( Block N " << current_block << " ) is diagonalized." << std::endl;
    }
 computeGroundEnergy();
}

void Hamiltonian::dump()
{
  BlockNumber NumberOfBlocks = parts.size();
  for (BlockNumber current_block=0;current_block<S.NumberOfBlocks();current_block++)
    {
      parts[current_block]->dump();
    }
  std::cout << "Hamiltonian has been dumped." << std::endl;
}

void Hamiltonian::computeGroundEnergy()
{
  RealVectorType LEV(S.NumberOfBlocks());
  BlockNumber NumberOfBlocks = parts.size();
  for (BlockNumber current_block=0;current_block<NumberOfBlocks;current_block++)
  {
	  LEV(current_block)=parts[current_block]->getMinimumEigenvalue();
  }
  GroundEnergy=LEV.minCoeff();
}

RealType Hamiltonian::getGroundEnergy()
{
  return GroundEnergy;
};

void Hamiltonian::reduce(const RealType Cutoff)
{
    std::cout << "Performing EV cutoff at " << Cutoff << " level" << std::endl;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber current_block=0;current_block<NumberOfBlocks;current_block++)
    {
	parts[current_block]->reduce(GroundEnergy+Cutoff);
    }
};

void Hamiltonian::save(H5::CommonFG* FG) const
{
    H5::Group RootGroup(FG->createGroup("Hamiltonian"));

    HDF5Storage::saveReal(&RootGroup,"GroundEnergy",GroundEnergy);

    // Save parts
    BlockNumber NumberOfBlocks = parts.size();
    H5::Group PartsGroup = RootGroup.createGroup("parts");
    for(BlockNumber n = 0; n < NumberOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.createGroup(nStr.str().c_str());
	parts[n]->save(&PartGroup);
    }
}

void Hamiltonian::load(const H5::CommonFG* FG)
{
    H5::Group RootGroup(FG->openGroup("Hamiltonian"));  

    GroundEnergy = HDF5Storage::loadReal(&RootGroup,"GroundEnergy");

    // FIXME!
    //if(Status!=Prepared) prepare();

    H5::Group PartsGroup = RootGroup.openGroup("parts");
    BlockNumber NumberOfBlocks = parts.size();
    if(NumberOfBlocks != PartsGroup.getNumObjs())
	throw(H5::GroupIException("Hamiltonian::load()","Inconsistent number of stored parts."));

    for(BlockNumber n = 0; n < NumberOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.openGroup(nStr.str().c_str());
	parts[n]->load(&PartGroup);
    }
    Status = Computed;
}

