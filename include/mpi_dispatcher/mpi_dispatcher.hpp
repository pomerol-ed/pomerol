//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/mpi_dispatcher/mpi_dispatcher.hpp
** \brief Implementation of a master-slave computation using unblocked MPI Communication that allows to utilize
** master as a computation node
*/
#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP
#define POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP

#include <mpi.h>

#include <cstddef>
#include <map>
#include <stack>
#include <vector>

namespace pMPI {

enum WorkerTag : int { Pending, Work, Finish }; // tags for MPI communication
using JobId = int;
using WorkerId = int;

struct MPIWorker {
    MPI_Comm Comm;
    WorkerId const id;
    WorkerId const boss;

    MPIWorker(MPI_Comm const& comm, WorkerId boss);
    void receive_order();
    void report_job_done();
    bool is_finished();
    bool is_working();

    JobId current_job() const { return current_job_; };

protected:
    JobId current_job_;
    MPI_Request req;
    WorkerTag Status;
};

struct MPIMaster {
    MPI_Comm Comm;
    int Ntasks, Nprocs;

    std::stack<JobId> JobStack;
    std::stack<WorkerId> WorkerStack;

    std::map<JobId, WorkerId> DispatchMap;
    std::vector<JobId> task_numbers;

    std::vector<WorkerId> worker_pool;
    std::map<std::size_t, WorkerId> WorkerIndices;

    std::vector<MPI_Request> wait_statuses;
    std::vector<bool> workers_finish;

    MPIMaster(MPI_Comm const& comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers);
    MPIMaster(MPI_Comm const& comm, std::vector<JobId> task_numbers, bool include_boss = true);
    MPIMaster(MPI_Comm const& comm, std::size_t ntasks, bool include_boss = true);

    void order_worker(WorkerId worker_id, JobId job);
    void order();
    void check_workers();
    bool is_finished() const;

private:
    void fill_stack_();
};

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP
