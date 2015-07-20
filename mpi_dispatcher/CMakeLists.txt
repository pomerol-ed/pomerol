set (mpi_dispatcher_src 
    mpi_dispatcher.hpp mpi_dispatcher.cpp 
    )

add_library(mpi_dispatcher OBJECT ${mpi_dispatcher_src})

if (Testing)
    add_subdirectory(test)
endif(Testing)

