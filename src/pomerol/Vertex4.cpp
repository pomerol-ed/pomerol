//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/Vertex4.cpp
/// \brief Irreducible two-particle vertex in the Matsubara representation.
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/Vertex4.hpp"

namespace Pomerol {

Vertex4::Vertex4(TwoParticleGF const& Chi,
                 GreensFunction const& G13,
                 GreensFunction const& G24,
                 GreensFunction const& G14,
                 GreensFunction const& G23)
    : Thermal(Chi.beta), ComputableObject(), Chi(Chi), G13(G13), G24(G24), G14(G14), G23(G23), Storage(*this) {}

void Vertex4::compute(long NumberOfMatsubaras) {
    if(getStatus() >= Computed)
        return;
    Storage.fill(NumberOfMatsubaras);
    setStatus(Computed);
}

ComplexType Vertex4::value(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const {
    ComplexType Value = Chi(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta * G13(MatsubaraNumber1) * G24(MatsubaraNumber2);
    if(MatsubaraNumber2 == MatsubaraNumber3)
        Value -= beta * G14(MatsubaraNumber1) * G23(MatsubaraNumber2);

    return Value;
}

ComplexType Vertex4::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const {
    return Storage(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);
}

} // namespace Pomerol
