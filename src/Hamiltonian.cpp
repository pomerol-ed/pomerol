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
#include "MPIDispatcher.h"

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#endif

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

void Hamiltonian::diagonalize(const boost::mpi::communicator & comm)
{
    if (Status >= Diagonalized) return;
    BlockNumber NumberOfBlocks = parts.size();
    int rank = comm.rank();
    comm.barrier();
    int comm_size = comm.size(); 
    if (rank==0) { INFO("Calculating using " << comm_size << " procs."); };

    size_t ROOT = 0;
    std::unique_ptr<MPI::MPIMaster> disp;

    if (comm.rank() == ROOT) { 
        // prepare one Master on a root process for distributing parts.size() jobs
        std::vector<MPI::JobId> job_order(parts.size());
        for (size_t i=0; i<job_order.size(); i++) job_order[i] = i;
        //for (auto x : job_order) std::cout << x << " " << std::flush; std::cout << std::endl;
        std::sort(job_order.begin(), job_order.end(), [&](MPI::JobId l, MPI::JobId r){return parts[l]->getSize() >= parts[r]->getSize();});
        //for (auto x : job_order) std::cout << x << " " << std::flush; std::cout << std::endl;
        disp.reset(new MPI::MPIMaster(comm,parts.size(),job_order,true)); 
        disp->order();
    };
    comm.barrier();

    // Start calculating data
    for (MPI::MPIWorker worker(comm,ROOT);!worker.is_finished();) {
        if (rank == ROOT) disp->order(); 
        worker.receive_order(); 
        if (worker.is_working()) { // for a specific worker
            auto p = worker.current_job;
            std::cout << "["<<p+1<<"/"<<parts.size()<< "] Proc " << comm.rank() << " : " << std::flush;
            parts[p]->diagonalize(); 
            worker.report_job_done(); 
	        INFO("Hpart " << p << " (" << S.getQuantumNumbers(BlockNumber(p)) << ") is diagonalized.");
        };
        if (rank == ROOT) disp->check_workers(); // check if there are free workers 
    };
    // at this moment all communication is finished
    // Now spread the information, who did what.
    comm.barrier();
    std::map<MPI::JobId, MPI::WorkerId> job_map;
    if (rank == ROOT) { 
        job_map = disp -> DispatchMap; 
        std::vector<MPI::JobId> jobs(job_map.size());
        std::vector<MPI::WorkerId> workers(job_map.size());
        std::transform(job_map.begin(), job_map.end(), jobs.begin(), [](std::pair<MPI::JobId, MPI::WorkerId> x){return x.first;});
        boost::mpi::broadcast(comm, jobs, ROOT);
        std::transform(job_map.begin(), job_map.end(), workers.begin(), [](std::pair<MPI::JobId, MPI::WorkerId> x){return x.second;});
        boost::mpi::broadcast(comm, workers, ROOT);
        disp.release(); 
    } 
    else {
        std::vector<MPI::JobId> jobs(parts.size());
        boost::mpi::broadcast(comm, jobs, ROOT);
        std::vector<MPI::WorkerId> workers(parts.size());
        boost::mpi::broadcast(comm, workers, ROOT);
        for (size_t i=0; i<jobs.size(); i++) job_map[jobs[i]] = workers[i]; 
    }

    //DEBUG for (auto x:job_map) std::cout << rank << ":" << x.first << "->" << x.second << std::endl;

    // Start distributing data
    comm.barrier();
    for (size_t p = 0; p<parts.size(); p++) {
            if (rank == job_map[p]){
                if (parts[p]->Status != HamiltonianPart::Diagonalized) { 
                    ERROR ("Worker" << rank << " didn't calculate part" << p); 
                    throw (std::logic_error("Worker didn't calculate this part."));
                    };
                boost::mpi::broadcast(comm, parts[p]->H.data(), parts[p]->H.rows()*parts[p]->H.cols(), rank);
                boost::mpi::broadcast(comm, parts[p]->Eigenvalues.data(), parts[p]->H.rows(), rank);
                }
            else {
                parts[p]->Eigenvalues.resize(parts[p]->H.rows());
                boost::mpi::broadcast(comm, parts[p]->H.data(), parts[p]->H.rows()*parts[p]->H.cols(), job_map[p]);
                boost::mpi::broadcast(comm, parts[p]->Eigenvalues.data(), parts[p]->H.rows(), job_map[p]);
                parts[p]->Status = HamiltonianPart::Diagonalized;
                 };
            };
/*
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++)
    {
	    parts[CurrentBlock]->diagonalize();
	    INFO("Hpart " << CurrentBlock << " (" << S.getQuantumNumbers(CurrentBlock) << ") is diagonalized.");
    }
*/
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

RealVectorType Hamiltonian::getEigenValues() const
{
    RealVectorType out(S.getNumberOfStates());
    size_t i=0;
    for (BlockNumber CurrentBlock=0; CurrentBlock<S.NumberOfBlocks(); CurrentBlock++) {
        const RealVectorType& tmp = parts[CurrentBlock]->getEigenValues();
        std::copy(tmp.data(), tmp.data() + tmp.size(), out.data()+i);
        i+=tmp.size(); 
        }
    return out;
}


RealType Hamiltonian::getGroundEnergy() const
{
    return GroundEnergy;
}

#ifdef ENABLE_SAVE_PLAINTEXT
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
#endif

} // end of namespace Pomerol
