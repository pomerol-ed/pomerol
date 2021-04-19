cmake_minimum_required(VERSION 3.1.0)

project(pomerol CXX)
message("\n${PROJECT_NAME} configuration\n")
string(TOUPPER ${PROJECT_NAME} PROJECT_CAPNAME)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(POMEROL_VERSION 2.0)
include(CommonDefs)

no_source_builds()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
fix_rpath()
set_linking_type()
compiler_quirks()

# Enable/Disable and find OpenMP
option(POMEROL_USE_OPENMP "Use OpenMP" TRUE)
if (POMEROL_USE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP")
        set(POMEROL_USE_OPENMP TRUE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        link_libraries(${OpenMP_LIBRARIES})
    else(OpenMP_CXX_FOUND)
        set(POMEROL_USE_OPENMP FALSE)
        message(STATUS "OpenMP not found - disabling")
    endif(OpenMP_CXX_FOUND)
else(POMEROL_USE_OPENMP)
    message(STATUS "OpenMP disabled")
endif(POMEROL_USE_OPENMP)

configure_file("${PROJECT_SOURCE_DIR}/include/pomerol/first_include.h.in"
               "${PROJECT_BINARY_DIR}/include/pomerol/first_include.h")
install(FILES "${PROJECT_BINARY_DIR}/include/pomerol/first_include.h"
        DESTINATION include/pomerol)

# TODO
set (pomerol_sources
    mpi_dispatcher/mpi_dispatcher
    pomerol/Misc
#    pomerol/Lattice
#    pomerol/LatticePresets
    pomerol/Index
#    pomerol/OperatorPresets
    pomerol/StatesClassification
    pomerol/HamiltonianPart
    pomerol/Hamiltonian
    pomerol/Thermal
    pomerol/MonomialOperatorPart
    pomerol/MonomialOperator
    pomerol/FieldOperatorContainer
    pomerol/DensityMatrixPart
    pomerol/DensityMatrix
    pomerol/GreensFunctionPart
    pomerol/GreensFunction
    pomerol/GFContainer
    pomerol/TwoParticleGFPart
    pomerol/TwoParticleGF
    pomerol/TwoParticleGFContainer
#    pomerol/Vertex4
#    pomerol/SusceptibilityPart
#    pomerol/Susceptibility
#    pomerol/EnsembleAverage
)

add_this_package(${pomerol_sources})

# Dependencies
add_eigen3()
add_libcommute()
add_mpi()
add_boost()
add_testing()

# Build executables
option(Progs "Build executables" OFF)
if (Progs)
    add_subdirectory(prog)
endif (Progs)

set(DOXYFILE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/src/mpi_dispatcher")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/mpi_dispatcher")
gen_documentation()

# Install pkg-config file
configure_file(pomerol.pc.in pomerol.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerol.pc
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

# Configure cmake config
configure_file(pomerolConfig.cmake.in pomerolConfig.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerolConfig.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pomerol)

# Configure lmod file
configure_file(pomerol.lmod.in pomerol.lmod @ONLY)
