#ifndef POMEROL_TEST_UTILITY_H
#define POMEROL_TEST_UTILITY_H

#include <mpi_dispatcher/misc.hpp>

#include <iostream>
#include <string>

void print_section(const std::string& str)
{
    if(!pMPI::rank(MPI_COMM_WORLD)) {
        std::cout << std::string(str.size(), '=') << std::endl;
        std::cout << str << std::endl;
        std::cout << std::string(str.size(), '=') << std::endl;
    }
}

#endif // #ifndef POMEROL_TEST_UTILITY_H
