set (mpi_dispatcher_src 
    mpi_dispatcher.hpp mpi_dispatcher.cpp 
    )

add_library(mpi_dispatcher OBJECT ${mpi_dispatcher_src})
set_target_properties(mpi_dispatcher PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(DIRECTORY . DESTINATION include/mpi_dispatcher 
          FILES_MATCHING PATTERN "*.hpp"
          PATTERN "test" EXCLUDE)

if (Testing)
    add_subdirectory(test)
endif(Testing)

