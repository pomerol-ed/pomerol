//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


#include "FieldOperator.h"

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include "mpi_dispatcher/mpi_skel.hpp"

namespace Pomerol{
FieldOperator::FieldOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    ComputableObject(), IndexInfo(IndexInfo), S(S), H(H), Index(Index)
{}

CreationOperator::CreationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,S,H,Index)
{
    O = new Pomerol::OperatorPresets::Cdag(Index);
}

AnnihilationOperator::AnnihilationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index) :
    FieldOperator(IndexInfo,S,H,Index)
{
    O = new Pomerol::OperatorPresets::C(Index);
}

FieldOperator::BlocksBimap const& FieldOperator::getBlockMapping() const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return LeftRightBlocks;
}
 
FieldOperatorPart& FieldOperator::getPartFromRightIndex(BlockNumber in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(in)->second];
}

FieldOperatorPart& FieldOperator::getPartFromRightIndex(const QuantumNumbers& in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromRight.find(S.getBlockNumber(in))->second];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(BlockNumber in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(in)->second];
}

FieldOperatorPart& FieldOperator::getPartFromLeftIndex(const QuantumNumbers& in) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    return *parts[mapPartsFromLeft.find(S.getBlockNumber(in))->second];
}

const std::vector<FieldOperatorPart*>& FieldOperator::getParts()
{
    return parts;
}

void FieldOperator::compute(const boost::mpi::communicator& comm)
{
    if (Status < Prepared) throw (exStatusMismatch());
    if (Status >= Computed) return;

    if (!comm.rank()) INFO_NONEWLINE("Computing " << *O << " in eigenbasis of the Hamiltonian: ");
/*

    pMPI::mpi_skel<pMPI::ComputeWrap<FieldOperatorPart>> skel;
    skel.parts.resize(parts.size());
    for (size_t i=0; i<parts.size(); i++) { skel.parts[i] = pMPI::ComputeWrap<FieldOperatorPart>(*parts[i]);};
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true); // actual running - very costly

    int rank = comm.rank();
    int comm_size = comm.size(); 

    // Start distributing data
    comm.barrier();
         
    for (size_t p = 0; p<parts.size(); p++) {

        int nelem = 0;
        if (rank == job_map[p]) nelem = parts[p]->elementsColMajor.nonZeros(); 
        boost::mpi::broadcast(comm, nelem, job_map[p]);
        auto data = parts[p]->elementsColMajor.data(); // Eigen::internal::CompressedStorage 
        auto& value_start = data.value(0);
        auto& index_start = data.index(0);
        if (comm.rank()) parts[p]->elementsColMajor.resize(nelem);
        exit(0);

        //boost::mpi::broadcast(comm, parts[p]->elementsColMajor.data(), nelem, job_map[p]);
        if (rank == job_map[p]) { 
                parts[p]->Status = FieldOperatorPart::Computed;
             };
        };
    comm.barrier();
*/
    size_t Size = parts.size();
    for (size_t BlockIn = 0; BlockIn < Size; BlockIn++){
        INFO_NONEWLINE( (int) ((1.0*BlockIn/Size) * 100 ) << "  " << std::flush);
        parts[BlockIn]->compute();
    };
    INFO("");
    Status = Computed;
}

ParticleIndex FieldOperator::getIndex(void) const
{
    return Index;
}

void CreationOperator::prepare(void)
{
    if (Status >= Prepared) return;
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0; RightIndex<S.NumberOfBlocks(); RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        //DEBUG(RightIndex << "->" << LeftIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new CreationOperatorPart(IndexInfo, S,
                                    H.getPart(RightIndex),H.getPart(LeftIndex),Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("CreationOperator_" << Index <<": " << Size << " parts will be computed");
    Status = Prepared;
}

void AnnihilationOperator::prepare()
{
    if (Status >= Prepared) return;
    size_t Size = parts.size();
    for (BlockNumber RightIndex=0;RightIndex<S.NumberOfBlocks();RightIndex++){
        BlockNumber LeftIndex = mapsTo(RightIndex);
        if (LeftIndex.isCorrect()){
            FieldOperatorPart *Part = new AnnihilationOperatorPart(IndexInfo, S,
                                    H.getPart(RightIndex),H.getPart(LeftIndex), Index);
            parts.push_back(Part);
            mapPartsFromRight[RightIndex]=Size;
            mapPartsFromLeft[LeftIndex]=Size;
            LeftRightBlocks.insert(BlockMapping(LeftIndex,RightIndex));
            Size++;
        }
    }
    INFO("AnnihilationOperator_" << Index <<": " << Size << " parts will be computed");
    Status = Prepared;
}

BlockNumber FieldOperator::getRightIndex(BlockNumber LeftIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    
    BlocksBimap::left_const_iterator it =  LeftRightBlocks.left.find(LeftIndex);
    return (it != LeftRightBlocks.left.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

BlockNumber FieldOperator::getLeftIndex(BlockNumber RightIndex) const
{
    if (Status < Prepared) { ERROR("FieldOperator is not prepared yet."); throw (exStatusMismatch()); }
    
    BlocksBimap::right_const_iterator it =  LeftRightBlocks.right.find(RightIndex);
    return (it != LeftRightBlocks.right.end()) ? it->second : ERROR_BLOCK_NUMBER;
}

BlockNumber FieldOperator::mapsTo(BlockNumber RightIndex) const
{
    bool found=false;
    std::map<FockState, MelemType> result;
    const std::vector<FockState> &states=S.getFockStates(RightIndex);
    for (std::vector<FockState>::const_iterator state_it=states.begin(); state_it!=states.end() && !found; state_it++) {
        result = O->actRight(*state_it);
        found = (result.size()>0);
        }
    return (found)?S.getBlockNumber(result.begin()->first):ERROR_BLOCK_NUMBER;
}

QuantumNumbers FieldOperator::mapsTo(const QuantumNumbers& in) const 
{
    BlockNumber out = this->mapsTo(S.getBlockNumber(in));
    if ( out == ERROR_BLOCK_NUMBER) throw (QuantumNumbers::exWrongNumbers());
    return S.getQuantumNumbers(out);
}

} // end of namespace Pomerol
