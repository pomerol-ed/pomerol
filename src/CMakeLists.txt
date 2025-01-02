#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(SOURCES
    mpi_dispatcher/mpi_dispatcher.cpp
    pomerol/Misc.cpp
    pomerol/LatticePresets.cpp
    pomerol/StatesClassification.cpp
    pomerol/HamiltonianPart.cpp
    pomerol/Hamiltonian.cpp
    pomerol/MonomialOperatorPart.cpp
    pomerol/MonomialOperator.cpp
    pomerol/FieldOperatorContainer.cpp
    pomerol/DensityMatrixPart.cpp
    pomerol/DensityMatrix.cpp
    pomerol/GreensFunctionPart.cpp
    pomerol/GreensFunction.cpp
    pomerol/GFContainer.cpp
    pomerol/TwoParticleGFPart.cpp
    pomerol/TwoParticleGF.cpp
    pomerol/TwoParticleGFContainer.cpp
    pomerol/Vertex4.cpp
    pomerol/SusceptibilityPart.cpp
    pomerol/Susceptibility.cpp
    pomerol/ThreePointSusceptibility.cpp
    pomerol/ThreePointSusceptibilityPart.cpp
    pomerol/ThreePointSusceptibilityContainer.cpp
    pomerol/EnsembleAverage.cpp
)

get_target_property(libcommute_INCLUDE_PATH libcommute::libcommute
                    INTERFACE_INCLUDE_DIRECTORIES)

add_library(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME}
                           PRIVATE
                           ${PROJECT_SOURCE_DIR}/include
                           ${PROJECT_BINARY_DIR}/include)
target_include_directories(${PROJECT_NAME}
                           SYSTEM PUBLIC
                           ${EIGEN3_INCLUDE_DIR}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${Boost_INCLUDE_DIRS}
                           ${libcommute_INCLUDE_PATH})
set_property(TARGET ${PROJECT_NAME} APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PRIVATE libcommute::libcommute)
target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES})
target_compile_options(${PROJECT_NAME} PUBLIC ${MPI_CXX_COMPILE_FLAGS}
                                              ${MPI_C_COMPILE_FLAGS})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(OpenMP_CXX_FOUND)
    target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
endif()

# Install the main library

install(TARGETS ${PROJECT_NAME}
        EXPORT "${PROJECT_NAME}Targets"
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
