//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/catch2/catch2-main.cpp
/// \brief User-provided main() for Catch2 that takes care of MPI initialization/finalization.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <mpi.h>

// A custom main() that takes care of MPI initialization/finalization
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run(argc, argv);
    MPI_Finalize();
    return result;
}
