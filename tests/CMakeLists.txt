include_directories (. ../src)
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


