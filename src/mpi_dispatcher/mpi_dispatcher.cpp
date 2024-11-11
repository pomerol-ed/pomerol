//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/mpi_dispatcher/mpi_dispatcher.cpp
/// \brief A master-worker parallelization scheme using non-blocking MPI communications (implementation).
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include "mpi_dispatcher/mpi_dispatcher.hpp"
#include "mpi_dispatcher/misc.hpp"

#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace pMPI {

//
// Worker
//

MPIWorker::MPIWorker(MPI_Comm const& comm, int boss)
    : Comm(comm), id(rank(comm)), boss(boss), current_job_(-1), req(MPI_REQUEST_NULL), Status(pMPI::Pending) {
    MPI_Irecv(&current_job_, 1, MPI_INT, boss, MPI_ANY_TAG, comm, &req);
}

bool MPIWorker::is_finished() {
    return Status == pMPI::Finish;
}

bool MPIWorker::is_working() {
    return Status == pMPI::Work;
}

void MPIWorker::receive_order() {
    if(Status != pMPI::Pending)
        return;

    MPI_Status st;
    int req_completed = 0;
    MPI_Test(&req, &req_completed, &st);

    if(req_completed) {
        Status = pMPI::WorkerTag(st.MPI_TAG);
        MPI_Irecv(&current_job_, 1, MPI_INT, boss, MPI_ANY_TAG, Comm, &req);
        if(is_finished()) {
            MPI_Cancel(&req);
        }
    }
}

void MPIWorker::report_job_done() {
    MPI_Send(nullptr, 0, MPI_INT, boss, pMPI::Pending, Comm);
    Status = pMPI::Pending;
}

//
// Master
//

void MPIMaster::fill_stack_() {
    for(int i = Ntasks - 1; i >= 0; --i) {
        JobStack.emplace(task_numbers[i]);
    }
    for(int p = Nprocs - 1; p >= 0; --p) {
        WorkerIndices[worker_pool[p]] = p;
        WorkerStack.emplace(worker_pool[p]);
    }
}

bool MPIMaster::is_finished() const {
    int NFinished = std::accumulate(workers_finish.begin(), workers_finish.end(), 0, std::plus<int>());
    return NFinished == Nprocs;
}

inline std::vector<WorkerId> _autorange_workers(MPI_Comm const& comm, bool include_boss) {
    int comm_size = pMPI::size(comm);
    int comm_rank = pMPI::rank(comm);

    std::vector<WorkerId> out;
    int Nprocs = comm_size - int(!include_boss);
    if(!Nprocs)
        throw std::runtime_error("No workers to evaluate");
    for(std::size_t p = 0; p < comm_size; ++p) {
        if(include_boss || comm_rank != p) {
            out.emplace_back(p);
        }
    }
    return out;
}

inline std::vector<JobId> _autorange_tasks(std::size_t ntasks) {
    std::vector<JobId> out(ntasks);
    std::iota(out.begin(), out.end(), 0);
    return out;
}

MPIMaster::MPIMaster(MPI_Comm const& comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers)
    : Comm(comm),
      Ntasks(static_cast<int>(task_numbers.size())),
      Nprocs(static_cast<int>(worker_pool.size())),
      task_numbers(std::move(task_numbers)),
      worker_pool(std::move(worker_pool)),
      wait_statuses(Nprocs, MPI_REQUEST_NULL),
      workers_finish(Nprocs, false) {
    fill_stack_();
}

MPIMaster::MPIMaster(MPI_Comm const& comm, std::size_t ntasks, bool include_boss)
    : MPIMaster(comm, _autorange_workers(comm, include_boss), _autorange_tasks(ntasks)) {}

MPIMaster::MPIMaster(MPI_Comm const& comm, std::vector<JobId> task_numbers, bool include_boss)
    : MPIMaster(comm, _autorange_workers(comm, include_boss), std::move(task_numbers)) {}

void MPIMaster::order_worker(WorkerId worker, JobId job) {
    MPI_Send(&job, 1, MPI_INT, worker, pMPI::Work, Comm);
    DispatchMap[job] = worker;
    MPI_Irecv(nullptr, 0, MPI_INT, worker, pMPI::Pending, Comm, &wait_statuses[WorkerIndices[worker]]);
}

void MPIMaster::order() {
    while(!WorkerStack.empty() && !JobStack.empty()) {
        WorkerId const& worker = WorkerStack.top();
        JobId const& job = JobStack.top();
        order_worker(worker, job);
        WorkerStack.pop();
        JobStack.pop();
    }
}

void MPIMaster::check_workers() {
    for(std::size_t i = 0; i < Nprocs; ++i) {
        if(wait_statuses[i] == MPI_REQUEST_NULL)
            continue;

        int req_completed = 0;
        MPI_Test(&wait_statuses[i], &req_completed, MPI_STATUS_IGNORE);
        if(req_completed) {
            WorkerStack.push(worker_pool[i]);
        }
    }
    if(JobStack.empty() && WorkerStack.size() >= Nprocs) {
        for(std::size_t i = 0; i < Nprocs; ++i) {
            if(!workers_finish[i]) {
                MPI_Send(nullptr, 0, MPI_INT, worker_pool[i], pMPI::Finish, Comm);
                // to prevent double sending of Finish command that could overlap with other communication
                workers_finish[i] = true;
            }
        }
    }
}

} // namespace pMPI
