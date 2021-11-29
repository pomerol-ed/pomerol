//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file include/mpi_dispatcher/mpi_skel.hpp
/// \brief Utilities for MPI-parallelized calculation of computable objects.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_SKEL_HPP
#define POMEROL_INCLUDE_MPI_DISPATCHER_MPI_SKEL_HPP

#include "misc.hpp"
#include "mpi_dispatcher.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <tuple>
#include <vector>

namespace pMPI {

/// \addtogroup MPI
///@{

/// \brief Wrapper around a computable object that calls the compute() method
/// of the wrapped object and carries information about the complexity of a call to that method.
/// \tparam PartType Type of the wrapped object.
template <typename PartType> struct ComputeWrap {
    /// Reference to the wrapped object.
    PartType& x;
    /// Complexity of a call to x.compute().
    int complexity;

    ComputeWrap() = default;
    /// Constructor.
    /// \param[in] x Reference to the wrapped object.
    /// \param[in] complexity Complexity of a call to x.compute().
    explicit ComputeWrap(PartType& x, int complexity = 1) : x(x), complexity(complexity) {}
    /// Call compute() of the wrapped object \ref x.
    void run() { x.compute(); }
};

/// \brief Wrapper around a computable object that calls the prepare() method
/// of the wrapped object and carries information about the complexity of a call to that method.
/// \tparam PartType Type of the wrapped object.
template <typename PartType> struct PrepareWrap {
    /// Reference to the wrapped object.
    PartType& x;
    /// Complexity of a call to x.prepare().
    int complexity;

    PrepareWrap() = default;
    /// Constructor.
    /// \param[in] x Reference to the wrapped object.
    /// \param[in] complexity Complexity of a call to x.prepare().
    explicit PrepareWrap(PartType& x, int complexity = 1) : x(x), complexity(complexity) {}
    /// Call prepare() of the wrapped object \ref x.
    void run() { x.prepare(); }
};

/// \brief This structure carries a list of wrappers and uses the mpi_dispatcher mechanism
/// to distribute the wrappers over MPI ranks and to call run() for all of them in parallel.
/// \tparam WrapType Type of the wrappers, one of \ref PrepareWrap and \ref ComputeWrap.
template <typename WrapType> struct mpi_skel {
    /// List of wrappers
    std::vector<WrapType> parts;
    /// Distribute the stored wrappers over MPI ranks according to their complexity
    /// and call run() for each of the wrappers.
    /// \param[in] Comm MPI communicator.
    /// \param[in] VerboseOutput Print extra information about the parallelization process.
    /// \return A mapping from job IDs to worker IDs assigned to perform the jobs.
    std::map<pMPI::JobId, pMPI::WorkerId> run(MPI_Comm const& Comm, bool VerboseOutput = true);
};

template <typename WrapType>
std::map<pMPI::JobId, pMPI::WorkerId> mpi_skel<WrapType>::run(MPI_Comm const& Comm, bool VerboseOutput) {
    int comm_rank = pMPI::rank(Comm);
    int comm_size = pMPI::size(Comm);
    int const root = 0;
    MPI_Barrier(Comm);

    if(comm_rank == root) {
        std::cout << "Calculating " << parts.size() << " jobs using " << comm_size << " procs." << std::endl;
    }

    std::unique_ptr<pMPI::MPIMaster> disp;

    if(comm_rank == root) {
        // prepare one Master on a root process for distributing parts.size() jobs
        std::vector<pMPI::JobId> job_order(parts.size());
        std::iota(job_order.begin(), job_order.end(), 0);

        auto comp1 = [this](std::size_t l, std::size_t r) -> int {
            return (parts[l].complexity > parts[r].complexity);
        };
        std::sort(job_order.begin(), job_order.end(), comp1);
        disp.reset(new pMPI::MPIMaster(Comm, job_order, true));
    }

    MPI_Barrier(Comm);

    // Start calculating data
    for(pMPI::MPIWorker worker(Comm, root); !worker.is_finished();) {
        if(comm_rank == root)
            disp->order();
        worker.receive_order();
        if(worker.is_working()) { // for a specific worker
            JobId p = worker.current_job();
            if(VerboseOutput)
                std::cout << "[" << p + 1 << "/" << parts.size() << "] P" << comm_rank << " : part " << p << " ["
                          << parts[p].complexity << "] run;" << std::endl;
            parts[p].run();
            worker.report_job_done();
        }
        if(comm_rank == root)
            disp->check_workers(); // check if there are free workers
    }

    // at this moment all communication is finished
    MPI_Barrier(Comm);
    // Now spread the information, who did what.
    if(VerboseOutput && comm_rank == root)
        std::cout << "done." << std::endl;

    MPI_Barrier(Comm);
    std::map<pMPI::JobId, pMPI::WorkerId> job_map;
    if(comm_rank == root) {
        job_map = disp->DispatchMap;
        long n_jobs = job_map.size();
        std::vector<pMPI::JobId> jobs(n_jobs);
        std::vector<pMPI::WorkerId> workers(n_jobs);

        auto it = job_map.cbegin();
        for(int i = 0; i < n_jobs; ++i, ++it) {
            std::tie(jobs[i], workers[i]) = *it;
        }

        MPI_Bcast(&n_jobs, 1, MPI_LONG, root, Comm);
        MPI_Bcast(jobs.data(), n_jobs, MPI_INT, root, Comm);
        MPI_Bcast(workers.data(), n_jobs, MPI_INT, root, Comm);
    } else {
        long n_jobs;
        MPI_Bcast(&n_jobs, 1, MPI_LONG, root, Comm);
        std::vector<pMPI::JobId> jobs(n_jobs);
        MPI_Bcast(jobs.data(), n_jobs, MPI_INT, root, Comm);
        std::vector<pMPI::WorkerId> workers(n_jobs);
        MPI_Bcast(workers.data(), n_jobs, MPI_INT, root, Comm);
        for(std::size_t i = 0; i < n_jobs; ++i)
            job_map[jobs[i]] = workers[i];
    }
    return job_map;
}

///@}

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_SKEL_HPP
