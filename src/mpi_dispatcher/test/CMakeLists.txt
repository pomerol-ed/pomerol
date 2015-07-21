#include_directories (. ../src)
enable_testing()

# MPI tests. Do not link to pomerol
set (mpi_tests
mpi_basic_test
)
#mpi_dispatcher_test

foreach (test ${mpi_tests})
    set(test_src ${test}.cpp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_NANOSLEEP") #1. old ubuntu gcc -D_GLIBCXX_USE_NANOSLEEP
    add_executable(${test} ${test_src} $<TARGET_OBJECTS:mpi_dispatcher>)
    target_link_libraries(${test} ${LINK_ALL} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES}) 
endforeach(test)

foreach (np 1 2 4 8 16)
foreach (test ${mpi_tests})
    set (test_parameters ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} "./${test}" ${MPIEXEC_POSTFLAGS})
    add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
endforeach(test)
endforeach(np)



