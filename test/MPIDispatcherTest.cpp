#include <mpi_dispatcher/misc.hpp>
#include <mpi_dispatcher/mpi_dispatcher.hpp>

#include "catch2/catch-pomerol.hpp"

#include <chrono>
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

TEST_CASE("Test mpi_dispatcher", "[mpi_dispatcher]") {
    std::random_device rd;
    std::mt19937 gen(100000);
    size_t ROOT = 0;
    int rank = pMPI::rank(MPI_COMM_WORLD);

    SECTION("With MPIMaster") {
        std::uniform_real_distribution<double> dist(0, 0.1);

        MPIWorker worker(MPI_COMM_WORLD, ROOT);
        dumb_task_type dumb_task;
        int ntasks = 45;

        std::unique_ptr<MPIMaster> disp;

        if (rank == ROOT) {
            disp.reset(new MPIMaster(MPI_COMM_WORLD, ntasks, true));
            disp->order();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        while(!worker.is_finished()) {
            if (rank == ROOT) disp->order();
            worker.receive_order();
            if (worker.is_working()) {
                dumb_task(dist(gen), worker.current_job(), rank);
                worker.report_job_done();
            };
            if (rank == ROOT)
                INFO("--> stack size = " << disp->JobStack.size() <<
                     " --> worker stack size =" << disp->WorkerStack.size() <<
                     "\n");
            if (rank == ROOT)
                disp->check_workers();
        }
        if (rank == ROOT) disp.release();

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task.counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(dumb_task.counter == ntasks);
    }

    SECTION("Without MPIMaster") {
        std::uniform_real_distribution<double> dist(0, 0.001);

        dumb_task_type dumb_task;
        int ntasks = 45;

        if (rank == ROOT) {
            MPIMaster master(MPI_COMM_WORLD, ntasks, false);
            for (; !master.is_finished();) {
                master.order();
                master.check_workers();
            }
        }
        else {
            MPIWorker worker(MPI_COMM_WORLD, ROOT);
            for (; !worker.is_finished();) {
                worker.receive_order();
                if (worker.is_working()) {
                    dumb_task(dist(gen), worker.current_job(), rank);
                    worker.report_job_done();
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task.counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(dumb_task.counter == ntasks);
    }
}
