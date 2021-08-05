//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2014 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2014 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.

/** \file prog/anderson.cpp
** \brief Diagonalization of the Anderson impurity model (1 impurity coupled to a set of non-interacting bath sites)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

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
