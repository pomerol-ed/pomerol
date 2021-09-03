#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required (VERSION 3.1.0)
project (pomerol_example CXX)

# Disable build in source
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(pomerol)
message(STATUS "Pomerol includes : ${pomerol_INCLUDE_DIRS}")
message(STATUS "Pomerol libraries : ${pomerol_LIBRARIES}")

include_directories(${pomerol_INCLUDE_DIRS})

add_executable(example2site example2site.cpp)
target_link_libraries(example2site ${pomerol_LIBRARIES})
