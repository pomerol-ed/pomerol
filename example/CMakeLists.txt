include_directories (. ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    pomerol
)

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (examples 
example2site
)

foreach (example ${examples})
    set(example_src ${example}.cpp)
    add_executable(${example} ${example_src})
endforeach(example)
