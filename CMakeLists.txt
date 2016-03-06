#
# Pomerol CMake configuration file
#

cmake_minimum_required (VERSION 2.8.12)
# title and version of the project
project (pomerol CXX)
message ("\n${PROJECT_NAME} configuration\n")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(POMEROL_VERSION 1.0)
# check ./cmake/CommondDefs.cmake for all of the macro/functions below
include(CommonDefs)

# Disable build in source
no_source_builds()
# Print build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
# rpath fix
fix_rpath()

# C++11 (optional)
option(AutoSetCXX11 "Auto set c++11 flags" OFF)
if (AutoSetCXX11)
    set_cxx11()
endif (AutoSetCXX11)

#Use complex Matrix elements
option(use_complex "Use complex matrix elements" OFF)
if (use_complex)
    message(STATUS "Using complex matrix elements")
    add_definitions("-DPOMEROL_COMPLEX_MATRIX_ELEMENS")
endif (use_complex)

set (pomerol_sources 
    mpi_dispatcher/mpi_dispatcher
    pomerol/Misc 
    pomerol/Lattice 
    pomerol/LatticePresets 
    pomerol/Index 
    pomerol/IndexClassification 
    pomerol/Operator 
    pomerol/OperatorPresets 
    pomerol/IndexHamiltonian 
    pomerol/Symmetrizer 
    pomerol/StatesClassification 
    pomerol/HamiltonianPart 
    pomerol/Hamiltonian 
    pomerol/Thermal 
    pomerol/FieldOperatorPart 
    pomerol/FieldOperator 
    pomerol/FieldOperatorContainer 
    pomerol/DensityMatrixPart 
    pomerol/DensityMatrix 
    pomerol/GreensFunctionPart 
    pomerol/GreensFunction
    pomerol/GFContainer
    pomerol/TwoParticleGFPart 
    pomerol/TwoParticleGF
    pomerol/TwoParticleGFContainer 
    pomerol/Vertex4 
    pomerol/Vertex4Container 
)

configure_file(${CMAKE_SOURCE_DIR}/include/pomerol.h.in ${CMAKE_BINARY_DIR}/include/pomerol.h)
INSTALL(FILES ${CMAKE_BINARY_DIR}/include/pomerol.h DESTINATION include)

add_this_package(${pomerol_sources})

# Dependencies
add_eigen3()
add_mpi()
add_boost(mpi serialization)
add_testing()

# Build executables
option(Progs "Build executables" OFF)
if (Progs)
    message(STATUS "Building executables")
    add_progs(anderson)
endif (Progs)

set(DOXYFILE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/src/mpi_dispatcher")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/mpi_dispatcher")
gen_documentation()

# Install pkg-config file
configure_file(pomerol.pc.in pomerol.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerol.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

# Configure cmake config
configure_file(pomerolConfig.cmake.in pomerolConfig.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerolConfig.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pomerol)

# Configure lmod file
configure_file(pomerol.lmod.in pomerol.lmod @ONLY)

