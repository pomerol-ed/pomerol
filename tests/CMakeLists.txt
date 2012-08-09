include_directories (. ../jsoncpp/include ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    jsoncpp
    pomerol
)

enable_testing()

set (tests 
IndexTermTest
IndexPermutationTest 
CCdagOperatorTest 
NOperatorTest 
SzOperatorTest 
HamiltonianPartTest 
HamiltonianTest 
FieldOperatorPartTest
FieldOperatorTest )

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
endforeach(test)
