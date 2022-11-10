#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.11.0)
project(pomerol_example CXX)

# Detect libcommute & pomerol
find_package(libcommute CONFIG REQUIRED)
find_package(pomerol CONFIG REQUIRED)

# Print out a bit of information about pomerol
get_target_property(pomerol_INCLUDE_DIRS pomerol::pomerol
                    INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(pomerol_LIBRARIES pomerol::pomerol
                    INTERFACE_LINK_LIBRARIES)
message(STATUS "Pomerol includes : ${pomerol_INCLUDE_DIRS}")
message(STATUS "Pomerol libraries : ${pomerol_LIBRARIES}")

# Add an executable and link it to pomerol
add_executable(example2site example2site.cpp)
target_link_libraries(example2site pomerol::pomerol)
