#include <mpi_dispatcher/mpi_dispatcher.hpp>
#include <thread>
#include <random>
#include <iostream>

using namespace pMPI;

int dumb_task_counter;

void dumb_task(double seconds, int jobid, int rank) {
    std::cout << "[" << rank << "] running job " << jobid << " " << seconds << " seconds..." << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(int(seconds * 1000)));
    ++dumb_task_counter;
    std::cout << "done." << std::endl;
};

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    std::random_device rd;
    std::mt19937 gen(100000);
    std::uniform_real_distribution<double> dist(0, 0.001);
    size_t ROOT = 0;
    int rank = pMPI::rank(MPI_COMM_WORLD);

    try {
        int ntasks = 7;
        dumb_task_counter = 0;

        std::unique_ptr<MPIWorker> worker_ptr;

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
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task_counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(dumb_task_counter != ntasks) {
            std::cout << "ntasks = " << ntasks
                      << ", dumb_task_counter = "
                      << dumb_task_counter << std::endl;
            return EXIT_FAILURE;
        }
    } // end try
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    };

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return EXIT_SUCCESS;
}
