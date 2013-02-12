//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#include "Hamiltonian.h"

#include<boost/filesystem.hpp>

namespace Pomerol{

Hamiltonian::Hamiltonian(const IndexClassification &IndexInfo, const IndexHamiltonian& F, const StatesClassification &S):
    ComputableObject(Constructed), IndexInfo(IndexInfo), F(F), S(S)
{}

Hamiltonian::~Hamiltonian()
{
    for(std::vector<boost::shared_ptr<HamiltonianPart> >::iterator iter = parts.begin(); iter != parts.end(); iter++)
	    iter->reset();
}

void Hamiltonian::prepare()
{
    if (Status >= Prepared) return;
    BlockNumber NumberOfBlocks = S.NumberOfBlocks();
    parts.resize(NumberOfBlocks);


    for (BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; CurrentBlock++)
    {
	    parts[CurrentBlock] =  boost::make_shared<HamiltonianPart>(HamiltonianPart(IndexInfo,F, S, CurrentBlock));
	    parts[CurrentBlock]->prepare();
	    INFO_NONEWLINE("Hpart " << CurrentBlock << " (" << S.getQuantumNumbers(CurrentBlock) << ") is entered. ");
        INFO("Size = " << S.getBlockSize(CurrentBlock));
    }
    Status = Prepared;
}

void Hamiltonian::diagonalize()
{
    if (Status >= Diagonalized) return;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++)
    {
	    parts[CurrentBlock]->diagonalize();
	    INFO("Hpart " << CurrentBlock << " (" << S.getQuantumNumbers(CurrentBlock) << ") is diagonalized.");
    }
    computeGroundEnergy();
    Status = Diagonalized;
}

void Hamiltonian::reduce(const RealType Cutoff)
{
    std::cout << "Performing EV cutoff at " << Cutoff << " level" << std::endl;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++)
    {
	parts[CurrentBlock]->reduce(GroundEnergy+Cutoff);
    }
}

void Hamiltonian::computeGroundEnergy()
{
    RealVectorType LEV(S.NumberOfBlocks());
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++) {
	    LEV(CurrentBlock) = parts[CurrentBlock]->getMinimumEigenvalue();
        }
    GroundEnergy=LEV.minCoeff();
}

const HamiltonianPart& Hamiltonian::getPart(const QuantumNumbers &in) const
{
    return *parts[S.getBlockNumber(in)];
}

const HamiltonianPart& Hamiltonian::getPart(BlockNumber in) const
{
    return *parts[in];
}

RealType Hamiltonian::getEigenValue(QuantumState state) const
{
    InnerQuantumState InnerState = S.getInnerState(state);
    return getPart(S.getBlockNumber(state)).getEigenValue(InnerState);
}

RealType Hamiltonian::getGroundEnergy() const
{
    return GroundEnergy;
}

void Hamiltonian::save(H5::CommonFG* RootGroup) const
{
    H5::Group HRootGroup(RootGroup->createGroup("Hamiltonian"));

    HDF5Storage::saveReal(&HRootGroup,"GroundEnergy",GroundEnergy);

    // Save parts
    BlockNumber NumberOfBlocks = parts.size();
    H5::Group PartsGroup = HRootGroup.createGroup("parts");
    for(BlockNumber n = 0; n < NumberOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.createGroup(nStr.str().c_str());
	parts[n]->save(&PartGroup);
    }
}

bool Hamiltonian::savetxt(const boost::filesystem::path &path)
{
    BlockNumber NumberOfBlocks = parts.size();
    boost::filesystem::create_directory(path);
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++) {
        std::stringstream tmp;
        tmp << "part" << S.getQuantumNumbers(CurrentBlock);
        boost::filesystem::path out = path / boost::filesystem::path (tmp.str()); 
	    parts[CurrentBlock]->savetxt(out);
        }
    return true;
}

void Hamiltonian::load(const H5::CommonFG* RootGroup)
{
    H5::Group HRootGroup(RootGroup->openGroup("Hamiltonian"));  

    GroundEnergy = HDF5Storage::loadReal(&HRootGroup,"GroundEnergy");

    H5::Group PartsGroup = HRootGroup.openGroup("parts");
    BlockNumber NumberOfBlocks = parts.size();
    if(NumberOfBlocks != PartsGroup.getNumObjs())
	throw(H5::GroupIException("Hamiltonian::load()","Inconsistent number of stored parts."));

    for(BlockNumber n = 0; n < NumberOfBlocks; n++){
	std::stringstream nStr;
	nStr << n;
	H5::Group PartGroup = PartsGroup.openGroup(nStr.str().c_str());
	parts[n]->load(&PartGroup);
    }
}

} // end of namespace Pomerol
