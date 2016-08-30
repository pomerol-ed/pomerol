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
AndersonTest02
AndersonTest03
)
if(use_complex)
 list(APPEND tests AndersonComplexTest)
endif(use_complex)

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


