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
#include "Logger.h"

namespace Pomerol {

namespace MPI {

enum class WorkerTag : int { Pending, Work, Finish }; // tags for MPI communication
typedef int JobId;
typedef int ProcId;

struct MPIWorker 
{
    boost::mpi::communicator Comm;
    const ProcId id;
    const ProcId boss;

    WorkerTag Status;

    JobId current_job = -1;

    boost::mpi::request WorkReq, FinishReq;

    MPIWorker(const boost::mpi::communicator &comm, ProcId boss):
        Comm(comm),
        id(Comm.rank()),
        boss(boss),
        WorkReq(Comm.irecv(boss, int(WorkerTag::Work), current_job)),
        FinishReq(Comm.irecv(boss, int(WorkerTag::Finish))),
        Status(WorkerTag::Pending)
    {};

    JobId get_job() { WorkReq.wait(); Status = WorkerTag::Work; return current_job;};
    void report_job_done() { Comm.isend(boss,int(WorkerTag::Pending)); Status = WorkerTag::Pending; WorkReq = Comm.irecv(boss, int(WorkerTag::Work), current_job);};
    bool is_finished();
};

struct MPIMaster 
{
    boost::mpi::communicator Comm;
    size_t Ntasks, Nprocs;
    const ProcId id;

    std::stack<JobId> JobStack;
    std::stack<ProcId> WorkerStack;

    //std::map<JobId, ProcId> DispatchMap;
    std::vector<JobId> task_numbers;

    std::vector<ProcId> worker_pool;
    std::map<size_t, ProcId> WorkerIndices;

    std::vector<boost::mpi::request> wait_statuses;

    MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, bool include_boss = true );
    MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, std::vector<ProcId> worker_pool, std::vector<JobId> task_numbers );

    void order_worker(ProcId worker_id, JobId job);
    void order();
    void update();
};

} // end of namespace MPI

} // end of namespace Pomerol

# endif // endif :: #ifndef __INCLUDE_MPIDISPATCHER_H

