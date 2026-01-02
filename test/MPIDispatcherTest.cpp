//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2026 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/MPIDispatcherTest.cpp
/// \brief Test the MPIMaster/MPIWorker communication mechanism.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include <mpi_dispatcher/misc.hpp>
#include <mpi_dispatcher/mpi_dispatcher.hpp>

#include "catch2/catch-pomerol.hpp"

#include <chrono>
#include <cstddef>
#include <memory>
#include <random>
#include <thread>

using namespace pMPI;

struct dumb_task_type {
    int counter = 0;
    void operator()(double seconds, int jobid, int rank) {
        std::this_thread::sleep_for(std::chrono::milliseconds(int(seconds * 1000)));
        ++counter;
    }
};

// cppcheck-suppress syntaxError
TEST_CASE("Test mpi_dispatcher", "[mpi_dispatcher]") {
    std::mt19937 gen(100000);
    int comm_rank = pMPI::rank(MPI_COMM_WORLD);
    int const root = 0;

    // cppcheck-suppress syntaxError
    SECTION("With MPIMaster") {
        std::uniform_real_distribution<double> dist(0, 0.1);

        MPIWorker worker(MPI_COMM_WORLD, root);
        dumb_task_type dumb_task;
        int ntasks = 45;

        std::unique_ptr<MPIMaster> disp(comm_rank == root ? new MPIMaster(MPI_COMM_WORLD, ntasks, true) : nullptr);

        if(comm_rank == root) {
            disp->order();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        while(!worker.is_finished()) {
            if(comm_rank == root)
                disp->order();
            worker.receive_order();
            if(worker.is_working()) {
                dumb_task(dist(gen), worker.current_job(), comm_rank);
                worker.report_job_done();
            };
            if(comm_rank == root) {
                INFO("--> stack size = " << disp->JobStack.size()
                                         << " --> worker stack size =" << disp->WorkerStack.size() << "\n");
                disp->check_workers();
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task.counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(dumb_task.counter == ntasks);
    }

    SECTION("Without MPIMaster") {
        std::uniform_real_distribution<double> dist(0, 0.001);

        dumb_task_type dumb_task;
        int ntasks = 45;

        if(comm_rank == root) {
            MPIMaster master(MPI_COMM_WORLD, ntasks, false);
            for(; !master.is_finished();) {
                master.order();
                master.check_workers();
            }
        } else {
            MPIWorker worker(MPI_COMM_WORLD, root);
            for(; !worker.is_finished();) {
                worker.receive_order();
                if(worker.is_working()) {
                    dumb_task(dist(gen), worker.current_job(), comm_rank);
                    worker.report_job_done();
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task.counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(dumb_task.counter == ntasks);
    }
}
