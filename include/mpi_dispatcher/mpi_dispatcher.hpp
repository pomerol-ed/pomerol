//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/mpi_dispatcher/mpi_dispatcher.hpp
/// \brief A master-worker parallelization scheme using non-blocking MPI communications.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP
#define POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP

#include <mpi.h>

#include <cstddef>
#include <map>
#include <stack>
#include <vector>

namespace pMPI {

/// \addtogroup MPI
///@{

/// MPI message tags used in communications between the master and its workers.
enum WorkerTag : int {
    Pending, ///< A worker is waiting for a new job
    Work,    ///< Request a worker to do a job
    Finish   ///< Order a worker to shut down
};

/// ID of a job.
using JobId = int;
/// ID of a worker process.
using WorkerId = int;

/// Abstraction of an MPI worker process.
struct MPIWorker {
    /// MPI communicator.
    MPI_Comm Comm;
    /// Worker ID of this process.
    WorkerId const id;
    /// Rank of the master process.
    int const boss;

    /// Constructor.
    /// \param[in] Comm MPI communicator
    /// \param[in] Boss Rank of the master process.
    MPIWorker(MPI_Comm const& Comm, int Boss);
    /// Check if there is an outstanding order from the master.
    void receive_order();
    /// Notify the master about a job's completion.
    void report_job_done();
    /// Has this worker process finished execution.
    bool is_finished();
    /// Is a job being processed by this worker?
    bool is_working();

    /// Get the ID of the job currently assigned to this worker.
    JobId current_job() const { return current_job_; };

protected:
    /// ID of the job currently assigned to this worker.
    JobId current_job_;
    /// An MPI request handle used for non-blocking communications.
    MPI_Request req;
    /// Current state of this worker.
    WorkerTag Status;
};

/// Abstraction of an MPI master process.
struct MPIMaster {
    /// MPI communicator.
    MPI_Comm Comm;
    /// Total number of jobs.
    int Ntasks;
    /// Total number of worker processes.
    int Nprocs;

    /// Stack of the jobs yet to be assigned to a worker.
    std::stack<JobId> JobStack;
    /// Stack of currently pending workers.
    std::stack<WorkerId> WorkerStack;

    /// A mapping from job IDs to IDs of the workers assigned to perform the jobs.
    std::map<JobId, WorkerId> DispatchMap;
    /// A list of IDs of all jobs to be completed.
    std::vector<JobId> task_numbers;
    /// A list of IDs of all worker processes.
    std::vector<WorkerId> worker_pool;
    /// Worker IDs and their serial numbers from the [0; worker_pool.size()[ range.
    std::map<WorkerId, int> WorkerIndices;

    /// MPI request handles used to perform non-blocking communications with the workers.
    std::vector<MPI_Request> wait_statuses;
    /// Flags to mark workers that have been shut down.
    std::vector<bool> workers_finish;

    /// Constructor.
    /// \param[in] Comm MPI communicator.
    /// \param[in] worker_pool A list of IDs of all worker processes.
    /// \param[in] task_numbers A list of IDs of all jobs to be completed.
    MPIMaster(MPI_Comm const& Comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers);
    /// Constructor. This version generates a list of worker IDs automatically.
    /// \param[in] Comm MPI communicator.
    /// \param[in] task_numbers A list of IDs of all jobs to be completed.
    /// \param[in] include_boss If true, allocate one worker per one MPI rank in the communicator.
    ///                         Otherwise, skip the rank of the master process.
    MPIMaster(MPI_Comm const& Comm, std::vector<JobId> task_numbers, bool include_boss = true);
    /// Constructor. This version generates lists of job IDs and worker IDs automatically.
    /// \param[in] Comm MPI communicator.
    /// \param[in] ntasks The number of jobs to be completed.
    /// \param[in] include_boss If true, allocate one worker per one MPI rank in the communicator.
    ///                         Otherwise, skip the rank of the master process.
    MPIMaster(MPI_Comm const& Comm, std::size_t ntasks, bool include_boss = true);

    /// Request a worker process to perform a job.
    /// \param[in] worker_id ID of the worker process.
    /// \param[in] job ID of the job to be performed.
    void order_worker(WorkerId worker_id, JobId job);
    /// Request the next available worker to perform the next job from the job stack.
    void order();
    /// Check which workers have become available and which have been shut down.
    void check_workers();
    /// Have all the workers been shut down?
    bool is_finished() const;

private:
    // Implementation details
    void fill_stack_();
};

///@}

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_HPP
