#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(DOXYFILE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/src/mpi_dispatcher")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/pomerol")
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${CMAKE_SOURCE_DIR}/include/mpi_dispatcher")

# Generate documentation
set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src")
set(DOXYFILE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/include")
set(DOXYFILE_IN "${PROJECT_SOURCE_DIR}/doc/Doxyfile.in")
set(DOXYFILE_NAME ${PROJECT_NAME})
set(DOXYFILE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")

include(UseDoxygen)
