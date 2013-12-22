#include <boost/mpi.hpp>

#include "MPIDispatcher.h"
#include <thread>
#include <random>

using namespace Pomerol;
using namespace Pomerol::MPI;

void dumb_task(double seconds) {
    std::cout << "running " << seconds << " seconds..." << std::flush;
//    int exec_time = 
    std::this_thread::sleep_for(std::chrono::milliseconds(int(seconds*1000)));
    std::cout << "done." << std::endl;
};

int main(int argc, char* argv[])
{
    Log.setDebugging(true);
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(world.rank()*24);
    std::uniform_real_distribution<double> dist(0,0.5);

    //std::cout << "Proc " << world.rank() << std::flush;


    if (world.rank() == 0) { // test sending to the same process
        INFO("Test : isend to the same proc");
        int i=4;
        int j = 5;
        INFO("i = " << i << ";" << "j = " << j);
        world.isend(0,0, i);
        world.isend(0,0, j);
        boost::mpi::request req[2];
        req[0] = world.irecv(0,0,j);
        req[1] = world.irecv(0,0,i);
        if (!world.rank()) boost::mpi::wait_all(req, req+2);
        INFO("i = " << i << ";" << "j = " << j);
        if (j!=4 || i!=5) return EXIT_FAILURE;
        INFO("Success : isend to the same proc");
    };

    try {

    int ntasks = 10;
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
        int p = worker.get_job();
        dumb_task(dist(gen));
        worker.report_job_done(); 
    
        if (rank == ROOT) disp->update(); 
    };
    if (rank == ROOT) disp.release();
    } // end try
    catch (std::exception &e){ERROR(e.what());};

    return 0;
} 
