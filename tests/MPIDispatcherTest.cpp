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
    
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    std::random_device rd;
    std::mt19937 gen(world.rank()*24);
    std::uniform_real_distribution<double> dist(0,0.1);
    size_t ROOT = 0;
    auto rank = world.rank();

    try {

    MPIWorker worker(world,ROOT);
    int ntasks = 15;

    std::unique_ptr<MPIMaster> disp;

    if (world.rank() == ROOT) { 
        disp.reset(new MPIMaster(world,ntasks,true));
        disp->order();
    };
    world.barrier();

    for (;!worker.is_finished();) {
        if (rank == ROOT) disp->order(); 
        worker.receive_order(); 
        if (worker.is_working()) {
            dumb_task(dist(gen));
            worker.report_job_done(); 
        };
        if (rank == ROOT) disp->check_workers(); 
    };
    if (rank == ROOT) disp.release();

    } // end try
    catch (std::exception &e){ERROR(e.what());return EXIT_FAILURE;};

    world.barrier();
    if (world.size() > 1) {
        int t = 10;
        DEBUG(world.rank() << " " << t);
        if (world.rank() == 0)  {t=45; world.isend(1,1,t);}
        else if (world.rank()==1) { 
            auto req = world.irecv(0,1,t); 
            auto msg = req.wait(); 
        };
        world.barrier();
        DEBUG(world.rank() << " " << t);
    };

// do it again
    try {

    int ntasks = 9;
    MPIWorker worker(world,ROOT);

    std::unique_ptr<MPIMaster> disp;

    if (world.rank() == ROOT) { 
        disp.reset(new MPIMaster(world,ntasks,true));
        disp->order();
    };
    world.barrier();

    for (;!worker.is_finished();) {
        if (rank == ROOT) disp->order(); 
        worker.receive_order(); 
        if (worker.is_working()) {
            dumb_task(dist(gen));
            worker.report_job_done(); 
        };
        if (rank == ROOT) disp->check_workers(); 
    };
    if (rank == ROOT) disp.release();
    } // end try
    catch (std::exception &e){ERROR(e.what());return EXIT_FAILURE;};

    return EXIT_SUCCESS;
} 
