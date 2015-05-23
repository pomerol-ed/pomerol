#include <boost/mpi.hpp>

#include "mpi_dispatcher/mpi_dispatcher.hpp"

using namespace pMPI;

int main(int argc, char* argv[])
{
    
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;


    if (world.rank() == 0) { // test sending to the same process
        std::cout << "Test : isend to the same proc" << std::endl;
        int i=4;
        int j = 5;
        std::cout << "i = " << i << ";" << "j = " << j << std::endl;
        world.isend(0,0, i);
        world.isend(0,0, j);
        int t1, t2;
        boost::mpi::request req[2];
        req[0] = world.irecv(0,0,t1);
        req[1] = world.irecv(0,0,t2);
        if (!world.rank()) boost::mpi::wait_all(req, req+2);
        i = t2; j = t1;
        std::cout << "i = " << i << ";" << "j = " << j << std::endl;
        if (!world.rank() && (j!=4 || i!=5)) return EXIT_FAILURE;
        std::cout << "Success : isend to the same proc" << std::endl;
    };

    return EXIT_SUCCESS;
} 
