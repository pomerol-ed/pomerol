include_directories (${CMAKE_SOURCE_DIR}/include/pomerol)
enable_testing()

# Pomerol tests
set (tests
OperatorTest
IndexPermutationTest
CCdagOperatorTest
NOperatorTest
SzOperatorTest
HamiltonianPartTest01
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
AndersonTest02
AndersonTest03
AndersonComplexTest
EnsembleAverage1siteTest
Susceptibility1siteTest
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

set(mpi_tests mpi_dispatcher_test mpi_dispatcher_test_nomaster broadcast)
foreach (test ${mpi_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test}
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        pomerol
    )

    foreach (np 2 4 8 16)
        set(test_parameters ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS})
        add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
    endforeach(np)
endforeach(test)
