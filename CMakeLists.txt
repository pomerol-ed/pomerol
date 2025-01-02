#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.11.0)

project(pomerol CXX)
set(POMEROL_VERSION 2.1)
set(POMEROL_DESCRIPTION "An exact diagonalization library aimed at \
solving condensed matter models of interacting fermions")
set(POMEROL_URL "https://pomerol-ed.github.io/pomerol/")

message(STATUS "Configuring ${PROJECT_NAME} ${POMEROL_VERSION}")

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Disable in-source builds
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled."
            "Please use a separate build directory.")
endif()

# Build type (Release by default)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# RPATH fix
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Library linking type
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Compiler quirks
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR
   "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        message(STATUS "Adding -stdlib=libc++ flag")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

# Enable/disable and find OpenMP
option(USE_OPENMP "Use OpenMP" ON)
if(USE_OPENMP)
    find_package(OpenMP)
else()
    message(STATUS "OpenMP disabled")
endif(USE_OPENMP)

#
# Dependencies
#

# Eigen 3
find_package(Eigen3 3.1 REQUIRED)

# libcommute
find_package(libcommute 0.7.2 CONFIG)

if(libcommute_FOUND)
    message(STATUS "Found libcommute version ${libcommute_VERSION}: "
                   "${libcommute_CONFIG}")
else(libcommute_FOUND)
    include(FetchContent)

    set(libcommute_GIT_REPOSITORY "https://github.com/krivenko/libcommute.git")
    set(libcommute_GIT_TAG "8f51813a7bb287c905bae6ba7c5de77d886ac9cb") # v0.7.2

    FetchContent_Declare(
        libcommute
        GIT_REPOSITORY ${libcommute_GIT_REPOSITORY}
        GIT_TAG        ${libcommute_GIT_TAG}
    )

    FetchContent_GetProperties(libcommute)
    if(NOT libcommute_POPULATED)
        message(STATUS
                "Fetching libcommute@${libcommute_GIT_TAG} "
                "from ${libcommute_GIT_REPOSITORY}"
        )

        FetchContent_Populate(libcommute)
        message(STATUS
                "libcommute has been fetched into ${libcommute_SOURCE_DIR}"
        )

        # Suppress compilation of libcommute's examples
        set(EXAMPLES OFF CACHE INTERNAL "" FORCE)
        # Add libcommute as a subproject
        add_subdirectory("${libcommute_SOURCE_DIR}" "${libcommute_BINARY_DIR}")
    endif(NOT libcommute_POPULATED)

    # Add an alias for the `libcommute` target in the subproject
    add_library(libcommute::libcommute ALIAS libcommute)
endif(libcommute_FOUND)

# MPI
find_package(MPI 3 REQUIRED)
message(STATUS "MPI includes: ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI C++ libs: ${MPI_CXX_LIBRARIES}")
message(STATUS "MPI flags: ${MPI_CXX_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")

#
# Workaround for MPI implementations that do not properly support
# MPI_CXX_* datatypes
#
option(Use_MPI_C_datatypes
    "Workaround: Use MPI_C_* datatypes instead of similar MPI_CXX_* datatypes"
    OFF)
mark_as_advanced(Use_MPI_C_datatypes)
if(Use_MPI_C_datatypes)
    add_definitions(-DPOMEROL_MPI_BOOL=MPI_C_BOOL
                    -DPOMEROL_MPI_DOUBLE_COMPLEX=MPI_C_DOUBLE_COMPLEX)
else(Use_MPI_C_datatypes)
    add_definitions(-DPOMEROL_MPI_BOOL=MPI_CXX_BOOL
                    -DPOMEROL_MPI_DOUBLE_COMPLEX=MPI_CXX_DOUBLE_COMPLEX)
endif(Use_MPI_C_datatypes)

# Boost
find_package(Boost 1.54.0 REQUIRED)
message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}" )

#
# Set up static analysis tools
#

option(StaticAnalysis "Run static analysis tools on C++ sources" OFF)
include(cmake/StaticAnalysis.cmake)

#
# Build parts of the project
#

# Main library
add_subdirectory(include)
add_subdirectory(src)

# Build executables
option(Progs "Build executables" OFF)
if(Progs)
    add_subdirectory(prog)
endif(Progs)

# Enable unit tests
option(Testing "Enable testing" ON)
if(Testing)
    enable_testing()
    add_subdirectory(test)
endif(Testing)

# Build documentation
option(Documentation "Build documentation" ON)
if(Documentation)
    add_subdirectory(doc)
endif()

#
# Configure and install CMake configs
#

include(CMakePackageConfigHelpers)

# Prepare version file
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${POMEROL_VERSION}
    COMPATIBILITY SameMajorVersion)

# Prepare config file
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})

# Install version & config files
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        DESTINATION lib/cmake/${PROJECT_NAME})

# Install targets file
install(EXPORT "${PROJECT_NAME}Targets" NAMESPACE ${PROJECT_NAME}::
        DESTINATION lib/cmake/${PROJECT_NAME})

#
# Configure and install misc files
#

# Install pkg-config file
configure_file(pomerol.pc.in pomerol.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerol.pc
        DESTINATION lib/pkgconfig)

# Configure lmod file
configure_file(pomerol.lmod.in pomerol.lmod @ONLY)
