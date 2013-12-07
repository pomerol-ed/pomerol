#include <boost/mpi.hpp>

int main()
{
    boost::mpi::environment MpiEnv;
    boost::mpi::communicator world;
    return 0;
} 
