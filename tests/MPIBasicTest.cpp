#include <boost/mpi.hpp>

#include "MPIDispatcher.h"

using namespace Pomerol;
using namespace Pomerol::MPI;

int main(int argc, char* argv[])
{
    
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;


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

    return EXIT_SUCCESS;
} 
