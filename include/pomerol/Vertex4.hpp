//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file include/pomerol/Vertex4.h
** \brief Irreducible two-particle vertex in the Matsubara representation.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#ifndef POMEROL_INCLUDE_VERTEX4_H
#define POMEROL_INCLUDE_VERTEX4_H

#include "ComputableObject.hpp"
#include "GreensFunction.hpp"
#include "MatsubaraContainers.hpp"
#include "Misc.hpp"
#include "Thermal.hpp"
#include "TwoParticleGF.hpp"

namespace Pomerol {

class Vertex4 : public Thermal, public ComputableObject {

    TwoParticleGF const& Chi4;
    GreensFunction const& G13;
    GreensFunction const& G24;
    GreensFunction const& G14;
    GreensFunction const& G23;

    /** Storage for precomputed values. */
    mutable MatsubaraContainer4<Vertex4> Storage;
    friend class MatsubaraContainer4<Vertex4>;

public:
    Vertex4(TwoParticleGF const& Chi4,
            GreensFunction const& G13,
            GreensFunction const& G24,
            GreensFunction const& G14,
            GreensFunction const& G23);

    void compute(long NumberOfMatsubaras = 0);

    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
    ComplexType value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_VERTEX4_H
