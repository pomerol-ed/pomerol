#include <boost/mpi.hpp>

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc,argv);
    boost::mpi::communicator world;
    return 0;
} 
