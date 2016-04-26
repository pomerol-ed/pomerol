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


/** \file src/TwoParticleGF.cpp
** \brief Two-particle Green's function in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#include "pomerol/TwoParticleGF.h"
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>

#include "mpi_dispatcher/mpi_skel.hpp"

namespace Pomerol{

TwoParticleGF::TwoParticleGF(const StatesClassification& S, const Hamiltonian& H,
                const AnnihilationOperator& C1, const AnnihilationOperator& C2, 
                const CreationOperator& CX3, const CreationOperator& CX4,
                const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(),
    S(S), H(H), C1(C1), C2(C2), CX3(CX3), CX4(CX4), DM(DM),
    parts(0), Vanishing(true),
    KroneckerSymbolTolerance (std::numeric_limits<RealType>::epsilon()), 
    ReduceResonanceTolerance (1e-8),
    CoefficientTolerance (1e-16), 
    ReduceInvocationThreshold (1e5)
{
}

TwoParticleGF::~TwoParticleGF()
{
      for(std::vector<TwoParticleGFPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

BlockNumber TwoParticleGF::getLeftIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber RightIndex) const
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getLeftIndex(RightIndex);
        case 1: return C2.getLeftIndex(RightIndex);
        case 2: return CX3.getLeftIndex(RightIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

BlockNumber TwoParticleGF::getRightIndex(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex) const
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getRightIndex(LeftIndex);
        case 1: return C2.getRightIndex(LeftIndex);
        case 2: return CX3.getRightIndex(LeftIndex);
        default: return ERROR_BLOCK_NUMBER;
    }
}

const FieldOperatorPart& TwoParticleGF::OperatorPartAtPosition(size_t PermutationNumber, size_t OperatorPosition, BlockNumber LeftIndex) const
{
    switch(permutations3[PermutationNumber].perm[OperatorPosition]){
        case 0: return C1.getPartFromLeftIndex(LeftIndex);
        case 1: return C2.getPartFromLeftIndex(LeftIndex);
        case 2: return CX3.getPartFromLeftIndex(LeftIndex);
        default: assert(0);
    }
    throw std::logic_error("TwoParticleGF : could not find operator part");
}

void TwoParticleGF::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of CX4.
    FieldOperator::BlocksBimap const& CX4NontrivialBlocks = CX4.getBlockMapping();
    for(FieldOperator::BlocksBimap::right_const_iterator outer_iter = CX4NontrivialBlocks.right.begin();
        outer_iter != CX4NontrivialBlocks.right.end(); outer_iter++){ // Iterate over the outermost index.
            for(size_t p=0; p<6; ++p){ // Choose a permutation
                  BlockNumber LeftIndices[4];
                  LeftIndices[0] = outer_iter->first;
                  LeftIndices[3] = outer_iter->second;
                  LeftIndices[2] = getLeftIndex(p,2,LeftIndices[3]);
                  LeftIndices[1] = getRightIndex(p,0,LeftIndices[0]);
                  // < LeftIndices[0] | O_1 | LeftIndices[1] >
                  // < LeftIndices[1] | O_2 | getRightIndex(p,1,LeftIndices[1]) >
                  // < LeftIndices[2]| O_3 | LeftIndices[3] >
                  // < LeftIndices[3] | CX4 | LeftIndices[0] >
                  // Select a relevant 'world stripe' (sequence of blocks).
                  if(getRightIndex(p,1,LeftIndices[1]) == LeftIndices[2] && LeftIndices[1].isCorrect() && LeftIndices[2].isCorrect()){
                      // DEBUG
                      /*DEBUG("new part: "  << S.getBlockInfo(LeftIndices[0]) << " " 
                                          << S.getBlockInfo(LeftIndices[1]) << " "
                                          << S.getBlockInfo(LeftIndices[2]) << " "
                                          << S.getBlockInfo(LeftIndices[3]) << " "
                      <<"BlockNumbers part: "  << LeftIndices[0] << " " << LeftIndices[1] << " " << LeftIndices[2] << " " << LeftIndices[3]);
                      */
                      parts.push_back(new TwoParticleGFPart(
                            OperatorPartAtPosition(p,0,LeftIndices[0]),
                            OperatorPartAtPosition(p,1,LeftIndices[1]),
                            OperatorPartAtPosition(p,2,LeftIndices[2]),
                            (CreationOperatorPart&)CX4.getPartFromLeftIndex(LeftIndices[3]),
                            H.getPart(LeftIndices[0]), H.getPart(LeftIndices[1]), H.getPart(LeftIndices[2]), H.getPart(LeftIndices[3]),
                            DM.getPart(LeftIndices[0]), DM.getPart(LeftIndices[1]), DM.getPart(LeftIndices[2]), DM.getPart(LeftIndices[3]),
                      permutations3[p]));

                      (*parts.rbegin())->KroneckerSymbolTolerance = KroneckerSymbolTolerance;
                      (*parts.rbegin())->ReduceResonanceTolerance = ReduceResonanceTolerance;
                      (*parts.rbegin())->CoefficientTolerance = CoefficientTolerance;
                      (*parts.rbegin())->ReduceInvocationThreshold = ReduceInvocationThreshold;
                      (*parts.rbegin())->MultiTermCoefficientTolerance = MultiTermCoefficientTolerance;
                      }
            }
    } 
    if ( parts.size() > 0 ) { 
        Vanishing = false;
        INFO("TwoParticleGF(" << getIndex(0) << getIndex(1) << getIndex(2) << getIndex(3) << "): " << parts.size() << " parts will be calculated");
        }
    Status = Prepared;
}

bool TwoParticleGF::isVanishing(void) const
{
    return Vanishing;
}


// An mpi adapter to 1) compute 2pgf terms; 2) convert them to a MatsubaraContainer; 3) purge terms
struct ComputeAndClearWrap
{
    int complexity;
    void run(){p->compute(); p->fillContainer(data_); if (clear_) p->clear();}; 
    ComputeAndClearWrap(TwoParticleGFPart::MatsubaraContainer &x, TwoParticleGFPart &p, bool clear, int complexity = 1):
        data_(x),p(&p),clear_(clear), complexity(complexity){};
protected:
    TwoParticleGFPart::MatsubaraContainer& data_;
    TwoParticleGFPart *p;
    bool clear_;
};



void TwoParticleGF::compute(bool clear, const boost::mpi::communicator & comm)
{

    TwoParticleGFPart::MatsubaraContainer data(DM.beta);

    if (Status < Prepared) throw (exStatusMismatch());
    if (Status >= Computed) return;
    if (!Vanishing) {
        // Create a "skeleton" class with pointers to part that can call a compute method
        pMPI::mpi_skel<ComputeAndClearWrap> skel;
        skel.parts.reserve(parts.size());
        for (size_t i=0; i<parts.size(); i++) { 
            skel.parts.push_back(ComputeAndClearWrap(data, *parts[i], clear, 1));
            };
        std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true); // actual running - very costly
        int rank = comm.rank();
        int comm_size = comm.size(); 

        // Start distributing data
        //DEBUG(comm.rank() << getIndex(0) << getIndex(1) << getIndex(2) << getIndex(3) << " Start distributing data");
        comm.barrier();
         
        for (size_t p = 0; p<parts.size(); p++) {
            boost::mpi::broadcast(comm, parts[p]->NonResonantTerms, job_map[p]);
            boost::mpi::broadcast(comm, parts[p]->ResonantTerms, job_map[p]);
            if (rank == job_map[p]) { 
                parts[p]->Status = TwoParticleGFPart::Computed;
                 };
            };
        comm.barrier();
    };
    Status = Computed;
}

// size_t TwoParticleGF::getNumResonantTerms() const
// {
//     size_t num = 0;
//     for(std::vector<TwoParticleGFPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
//         num += (*iter)->getNumResonantTerms();
//     return num;
// }
// 
// size_t TwoParticleGF::getNumNonResonantTerms() const
// {
//     size_t num = 0;
//     for(std::vector<TwoParticleGFPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
//         num += (*iter)->getNumNonResonantTerms();
//     return num;
// }

ParticleIndex TwoParticleGF::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C1.getIndex();
        case 1: return C2.getIndex();
        case 2: return CX3.getIndex();
        case 3: return CX4.getIndex();
        default: assert(0);
    }
    throw std::logic_error("TwoParticleGF : could not get operator index");
}

unsigned short TwoParticleGF::getPermutationNumber ( const Permutation3& in )
{
    for (unsigned short i=0; i<6; ++i) if (in == permutations3[i]) return i;
    ERROR("TwoParticleGF: Permutation " << in << " not found in all permutations3");
    return 0;
}

} // end of namespace Pomerol

