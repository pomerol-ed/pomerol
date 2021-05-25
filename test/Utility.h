#include <mpi_dispatcher/misc.hpp>

#include <iostream>
#include <string>

// Generalized 'square' function.
template<typename T> inline T sqr(T x) { return x*x; }

void print_section(const std::string& str)
{
    if(!pMPI::rank(MPI_COMM_WORLD)) {
        std::cout << std::string(str.size(), '=') << std::endl;
        std::cout << str << std::endl;
        std::cout << std::string(str.size(), '=') << std::endl;
    }
}
