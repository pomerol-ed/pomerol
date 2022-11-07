#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Build Catch2 main object file
add_library(catch2 OBJECT catch2/catch2-main.cpp)
target_include_directories(catch2 PRIVATE ${MPI_CXX_INCLUDE_PATH})
if(StaticAnalysis)
    set_target_properties(catch2 PROPERTIES CXX_CLANG_TIDY "" CXX_CPPCHECK "")
endif()

# Pomerol tests
set(tests
    HamiltonianTest
    HamiltonianBosonsTest
    GF1siteTest
    GF2siteTest
    GF4siteTest
    AndersonTest
    AndersonComplexTest
    Anderson2PGFTest
    Vertex4Test
    SusceptibilityTest
    3PSusc1siteTest
    3PSusc3siteTest
)

foreach(test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    set(test_parameters
        ${MPIEXEC_NUMPROC_FLAG} 1
        ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS})
    add_test(NAME ${test} COMMAND "${MPIEXEC}" ${test_parameters})
    target_link_libraries(${test} PRIVATE
                          ${PROJECT_NAME} ${MPI_CXX_LIBRARIES} catch2)
endforeach(test)

set(mpi_tests BroadcastTest MPIDispatcherTest)
foreach(test ${mpi_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} PRIVATE
                          ${PROJECT_NAME} ${MPI_CXX_LIBRARIES} catch2)

    foreach(np 2 4 8 16)
        set(test_parameters
            ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS})
        add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
    endforeach(np)
endforeach(test)
