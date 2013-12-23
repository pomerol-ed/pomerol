#include <boost/mpi.hpp>

#include "MPIDispatcher.h"
#include <thread>
#include <random>

using namespace Pomerol;
using namespace Pomerol::MPI;

void dumb_task(double seconds) {
    std::cout << "running " << seconds << " seconds..." << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(int(seconds*1000)));
    std::cout << "done." << std::endl;
};

int main(int argc, char* argv[])
{
    Log.setDebugging(true);
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    std::random_device rd;
    std::mt19937 gen(world.rank()*24);
    std::uniform_real_distribution<double> dist(0,0.25);

    try {

    int ntasks = 15;
    size_t ROOT = 0;
    auto rank = world.rank();

    std::unique_ptr<MPIMaster> disp;

    if (world.rank() == ROOT) { 
        disp.reset(new MPIMaster(world,ntasks,true));
        disp->order();
    };
    world.barrier();

    for (MPIWorker worker(world,ROOT);!worker.is_finished();) {
        if (rank == ROOT) disp->order(); 
        worker.receive_order(); 
        if (worker.is_working()) {
            dumb_task(dist(gen));
            worker.report_job_done(); 
        };
        if (rank == ROOT) disp->update(); 
    };
    if (rank == ROOT) disp.release();
    } // end try
    catch (std::exception &e){ERROR(e.what());return EXIT_FAILURE;};

    return EXIT_SUCCESS;
} 
