#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

if(StaticAnalysis)

    # clang-tidy
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.6)
        if(NOT CMAKE_CXX_CLANG_TIDY)
            find_program(CLANG_TIDY_EXECUTABLE
                         NAMES "clang-tidy" REQUIRED
                         DOC "Path to the clang-tidy executable")
            message(STATUS "Using clang-tidy: ${CLANG_TIDY_EXECUTABLE}")
            set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
        endif()
    else()
        message(WARNING "CMake >= 3.6 is required to run clang-tidy checks")
    endif()

    # cppcheck
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.10)
        if(NOT CMAKE_CXX_CPPCHECK)
            find_program(CPPCHECK_EXECUTABLE
                         NAMES "cppcheck" REQUIRED
                         DOC "Path to the cppcheck executable")
            message(STATUS "Using cppcheck: ${CPPCHECK_EXECUTABLE}")

            configure_file(.cppcheck.supp.in .cppcheck.supp @ONLY)

            set(CMAKE_CXX_CPPCHECK
                "${CPPCHECK_EXECUTABLE}"
                 "--enable=warning,style,performance,portability"
                 "--std=c++11"
                 "--template=gcc"
                 "--suppressions-list=${CMAKE_BINARY_DIR}/.cppcheck.supp"
                 "--inline-suppr"
                 "--verbose"
                 "--force"
                 "--quiet"
            )
        endif()
    else()
        message(WARNING "CMake >= 3.10 is required to run cppcheck")
    endif()

endif()
