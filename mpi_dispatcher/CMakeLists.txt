set (mpi_dispatcher_src 
    mpi_dispatcher.hpp mpi_dispatcher.cpp 
    )

add_library(mpi_dispatcher STATIC ${mpi_dispatcher_src})

if (Testing)
    add_subdirectory(test)
endif(Testing)

#add_subdirectory(prog)


