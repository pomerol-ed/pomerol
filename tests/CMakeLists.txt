include_directories (. ../src)
enable_testing()

# MPI tests. Do not link to pomerol
set (mpi_tests
MPIBasicTest
MPIDispatcherTest
)

foreach (test ${mpi_tests})
    set(test_src ${test}.cpp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_NANOSLEEP") #1. old ubuntu gcc -D_GLIBCXX_USE_NANOSLEEP
    add_executable(${test} ${test_src} ${CMAKE_SOURCE_DIR}/src/MPIDispatcher.cpp)
    target_link_libraries(${test} 
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
    )
endforeach(test)

foreach (np 1 2 4 8 16)
foreach (test ${mpi_tests})
    set (test_parameters -np ${np} "./${test}")
    add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
endforeach(test)
endforeach(np)


# Pomerol tests
set (tests 
OperatorTest
IndexPermutationTest 
CCdagOperatorTest 
NOperatorTest 
SzOperatorTest 
HamiltonianPartTest01 
#HamiltonianPartTest02 
#SingletTest
HamiltonianTest 
FieldOperatorPartTest
FieldOperatorTest 
GF1siteTest
GF2siteTest
AndersonTest
GF4siteTest
GFContainerTest
TwoParticleGFContainerTest
Vertex4Test
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    set (test_parameters -np 1 "./${test}")
    add_test(NAME ${test} COMMAND "${MPIEXEC}" ${test_parameters})
    target_link_libraries(${test} 
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        pomerol
    )
endforeach(test)


set (tests_np 
)

foreach (np 1 4 16)
foreach (test ${tests_np})
    set (test_parameters -np ${np} "./${test}")
    add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
endforeach(test)
endforeach(np)


