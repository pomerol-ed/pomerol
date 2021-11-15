//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/mpi_dispatcher/mpi_skel.hpp
** \brief Declares mpi_skel - a structure to simplify master-slave calculations
*/
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

template <typename PartType> struct ComputeWrap {
    PartType& x;
    int complexity;

    ComputeWrap() = default;
    explicit ComputeWrap(PartType& x, int complexity = 1) : x(x), complexity(complexity) {}
    void run() { x.compute(); }
};

template <typename PartType> struct PrepareWrap {
    PartType& x;
    int complexity;

    PrepareWrap() = default;
    explicit PrepareWrap(PartType& x, int complexity = 1) : x(x), complexity(complexity) {}
    void run() { x.prepare(); }
};

template <typename WrapType> struct mpi_skel {
    std::vector<WrapType> parts;
    std::map<pMPI::JobId, pMPI::WorkerId> run(MPI_Comm const& comm, bool VerboseOutput = true);
};

/// Master-slave task schedule, associates jobs with workers and executes mpi_dispatcher to run the simulation
template <typename WrapType>
std::map<pMPI::JobId, pMPI::WorkerId> mpi_skel<WrapType>::run(MPI_Comm const& comm, bool VerboseOutput) {
    int comm_rank = pMPI::rank(comm);
    int comm_size = pMPI::size(comm);
    int const root = 0;
    MPI_Barrier(comm);

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
        disp.reset(new pMPI::MPIMaster(comm, job_order, true));
    }

    MPI_Barrier(comm);

    // Start calculating data
    for(pMPI::MPIWorker worker(comm, root); !worker.is_finished();) {
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
    MPI_Barrier(comm);
    // Now spread the information, who did what.
    if(VerboseOutput && comm_rank == root)
        std::cout << "done." << std::endl;

    MPI_Barrier(comm);
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

        MPI_Bcast(&n_jobs, 1, MPI_LONG, root, comm);
        MPI_Bcast(jobs.data(), n_jobs, MPI_INT, root, comm);
        MPI_Bcast(workers.data(), n_jobs, MPI_INT, root, comm);
    } else {
        long n_jobs;
        MPI_Bcast(&n_jobs, 1, MPI_LONG, root, comm);
        std::vector<pMPI::JobId> jobs(n_jobs);
        MPI_Bcast(jobs.data(), n_jobs, MPI_INT, root, comm);
        std::vector<pMPI::WorkerId> workers(n_jobs);
        MPI_Bcast(workers.data(), n_jobs, MPI_INT, root, comm);
        for(std::size_t i = 0; i < n_jobs; ++i)
            job_map[jobs[i]] = workers[i];
    }
    return job_map;
}

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_SKEL_HPP
