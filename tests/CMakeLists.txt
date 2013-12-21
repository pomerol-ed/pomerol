include_directories (. ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    pomerol
)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (tests 
MPITest
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
endforeach(test)
