#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Install header files
install(DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.hpp")

# Prepare and install pomerol/Version.hpp
if(USE_OPENMP AND OPENMP_FOUND)
    set(POMEROL_USE_OPENMP ON)
endif()
configure_file("pomerol/Version.hpp.in" "pomerol/Version.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pomerol/Version.hpp"
        DESTINATION include/pomerol)
