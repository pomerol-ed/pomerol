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


/** \file src/DensityMatrix.h
** \brief A storage of the matrix elements of the hamiltonian in Fock basis, provides eigenvalues and eigenfunctions
** 
** \author Andrey Antipov(Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_MPIDISPATCHER_H
#define __INCLUDE_MPIDISPATCHER_H

#include <boost/mpi.hpp>
#include <stack>

#include "Misc.h"

namespace Pomerol {

namespace pMPI {

//enum class WorkerTag : int { Pending, Work, Finish }; // tags for MPI communication
enum WorkerTag { Pending, Work, Finish }; // tags for MPI communication
typedef int JobId;
typedef int WorkerId;

struct MPIWorker 
{
    boost::mpi::communicator Comm;
    const WorkerId id;
    const WorkerId boss;

    WorkerTag Status;
    boost::mpi::request WorkReq, FinishReq;
    JobId current_job;

    MPIWorker(const boost::mpi::communicator &comm, WorkerId boss);
    void receive_order();
    void report_job_done();
    bool is_finished();
    bool is_working();
};

struct MPIMaster 
{
    boost::mpi::communicator Comm;
    size_t Ntasks, Nprocs;
    const WorkerId id;

    std::stack<JobId> JobStack;
    std::stack<WorkerId> WorkerStack;

    std::map<JobId, WorkerId> DispatchMap;
    std::vector<JobId> task_numbers;

    std::vector<WorkerId> worker_pool;
    std::map<size_t, WorkerId> WorkerIndices;

    std::vector<boost::mpi::request> wait_statuses;
    std::vector<bool> workers_finish;

    void swap(MPIMaster &x);

    MPIMaster(const boost::mpi::communicator &comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers );
    MPIMaster(const boost::mpi::communicator &comm, std::vector<JobId> task_numbers, bool include_boss = true );
    MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, bool include_boss = true );

    void order_worker(WorkerId worker_id, JobId job);
    void order();
    void check_workers();
protected:
    void fill_stack_();
};

} // end of namespace MPI

} // end of namespace Pomerol

# endif // endif :: #ifndef __INCLUDE_MPIDISPATCHER_H

