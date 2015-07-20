set (mpi_dispatcher_src 
    mpi_dispatcher.hpp mpi_dispatcher.cpp 
    )

add_library(mpi_dispatcher STATIC ${mpi_dispatcher_src})
set_target_properties(mpi_dispatcher PROPERTIES POSITION_INDEPENDENT_CODE ON)
install(TARGETS mpi_dispatcher 
        EXPORT mpi_dispatcher
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        )
install(EXPORT mpi_dispatcher NAMESPACE pomerol:: DESTINATION share/pomerol/mpi_dispatcher)

if (Testing)
    add_subdirectory(test)
endif(Testing)

