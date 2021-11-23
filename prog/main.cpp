//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file prog/main.cpp
/// \brief A common main source file for executables that diagonalize various quantum models.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#ifdef POMEROL_ANDERSON
#include "anderson_model.hpp"
#endif
#ifdef POMEROL_HUBBARD2D
#include "hubbard2d_model.hpp"
#endif

int main(int argc, char* argv[]) {
#ifdef POMEROL_ANDERSON
    anderson_model m(argc, argv);
#endif
#ifdef POMEROL_HUBBARD2D
    hubbard2d_model m(argc, argv);
#endif

    m.compute();
    return 0;
}
