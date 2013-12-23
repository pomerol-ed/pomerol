include_directories (. ../src)
enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
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
Vertex4Test
#Vertex4ContainerTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    set (test_parameters -np 1 "./${test}")
    add_test(NAME ${test} COMMAND "${MPIEXEC}" ${test_parameters})
    target_link_libraries(${test} 
        ${Boost_LIBRARIES}
#        ${HDF5_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        pomerol
    )
endforeach(test)

set (mpi_tests
MPIBasicTest
MPIDispatcherTest
)

foreach (test ${mpi_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} 
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        pomerol
    )
endforeach(test)

foreach (np 1 2 4 8 16)
foreach (test ${mpi_tests})
    set (test_parameters -np ${np} "./${test}")
    add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
endforeach(test)
endforeach(np)

