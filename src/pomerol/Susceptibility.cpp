//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/Susceptibility.cpp
/// \brief Dynamical susceptibility in the Matsubara representation (implementation).
/// \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/Susceptibility.hpp"

namespace Pomerol {

Susceptibility::Susceptibility(StatesClassification const& S,
                               Hamiltonian const& H,
                               MonomialOperator const& A,
                               MonomialOperator const& B,
                               DensityMatrix const& DM)
    : Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), B(B), DM(DM) {}

Susceptibility::Susceptibility(Susceptibility const& Chi)
    : Thermal(Chi.beta),
      ComputableObject(Chi),
      S(Chi.S),
      H(Chi.H),
      A(Chi.A),
      B(Chi.B),
      DM(Chi.DM),
      Vanishing(Chi.Vanishing),
      parts(Chi.parts),
      SubtractDisconnected(Chi.SubtractDisconnected),
      ave_A(Chi.ave_A),
      ave_B(Chi.ave_B) {}

void Susceptibility::prepare() {
    if(getStatus() >= Prepared)
        return;

    // Find out non-trivial blocks of A and B.
    auto const& ANontrivialBlocks = A.getBlockMapping();
    auto const& BNontrivialBlocks = B.getBlockMapping();

    auto Aiter = ANontrivialBlocks.left.begin();
    auto Biter = BNontrivialBlocks.right.begin();

    while(Aiter != ANontrivialBlocks.left.end() && Biter != BNontrivialBlocks.right.end()) {
        // <Aleft|A|Aright><Bleft|B|Bright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;
        BlockNumber Bleft = Biter->second;
        BlockNumber Bright = Biter->first;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(Aleft == Bright && Aright == Bleft) {
            // check if retained blocks are included. If not, do not push.
            if(DM.isRetained(Aleft) || DM.isRetained(Aright))
                parts.emplace_back(const_cast<MonomialOperatorPart&>(A.getPartFromLeftIndex(Aleft)),
                                   const_cast<MonomialOperatorPart&>(B.getPartFromRightIndex(Bright)),
                                   H.getPart(Aright),
                                   H.getPart(Aleft),
                                   DM.getPart(Aright),
                                   DM.getPart(Aleft));
        }

        unsigned long AleftInt = Aleft;
        unsigned long BrightInt = Bright;

        if(AleftInt <= BrightInt)
            Aiter++;
        if(AleftInt >= BrightInt)
            Biter++;
    }

    if(!parts.empty())
        Vanishing = false;

    setStatus(Prepared);
}

void Susceptibility::compute() {
    if(getStatus() >= Computed)
        return;
    if(getStatus() < Prepared)
        prepare();

    if(getStatus() < Computed) {
        for(auto& p : parts)
            p.compute();
    }
    setStatus(Computed);
}

void Susceptibility::subtractDisconnected() {
    EnsembleAverage EA_A(A, DM);
    EnsembleAverage EA_B(B, DM);
    subtractDisconnected(EA_A, EA_B);
}

void Susceptibility::subtractDisconnected(ComplexType ave_A, ComplexType ave_B) {
    SubtractDisconnected = true;
    this->ave_A = ave_A;
    this->ave_B = ave_B;
}

void Susceptibility::subtractDisconnected(EnsembleAverage& EA_A, EnsembleAverage& EA_B) {
    EA_A.compute();
    EA_B.compute();
    subtractDisconnected(EA_A(), EA_B());
}

} // namespace Pomerol
