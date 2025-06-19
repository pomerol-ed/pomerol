//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/GF1siteTest.cpp
/// \brief Single-particle anomalous Green's functions of a single Hubbard atom.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GFContainer.hpp>
#include <pomerol/GreensFunction.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/StatesClassification.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>
#include <set>

using namespace Pomerol;

// cppcheck-suppress syntaxError
TEST_CASE("Anomalous Green's function of a Hubbard atom", "[GF1siteAn]") {
    RealType U = 1.0;
    RealType mu = 0.4;
    RealType Delta = 0.1;
    RealType beta = 10.0;

    RealType e1 = -mu;
    RealType e2 = -2 * mu + U;
    RealType ep = std::sqrt(e2 * e2 + 4 * Delta * Delta);
    RealType esc1 = 0.5 * (e2 - ep);
    RealType esc2 = 0.5 * (e2 + ep);

    RealType a = esc2 / std::sqrt(Delta * Delta + esc2 * esc2);
    RealType b = esc1 / std::sqrt(Delta * Delta + esc1 * esc1);

    RealType w1 = std::exp(-beta * e1);
    RealType wsc1 = std::exp(-beta * esc1);
    RealType wsc2 = std::exp(-beta * esc2);
    RealType Z = 2 * w1 + wsc1 + wsc2;
    w1 /= Z;
    wsc1 /= Z;
    wsc2 /= Z;

    // Reference normal Green's function
    auto G_ref = [&](int n) {
        RealType omega = M_PI * (2 * n + 1) / beta;
        return (w1 + wsc1) * (a * a) / (I * omega - (e1 - esc1)) +
               (w1 + wsc2) * (1 - a * a) / (I * omega - (e1 - esc2)) +
               (wsc1 + w1) * (b * b) / (I * omega - (esc1 - e1)) +
               (wsc2 + w1) * (1 - b * b) / (I * omega - (esc2 - e1));
    };
    // Reference anomalous Green's function
    auto F_ref = [&](int n) {
        RealType omega = M_PI * (2 * n + 1) / beta;
        return (Delta / ep) * (-(w1 + wsc1) / (I * omega - (e1 - esc1)) + (w1 + wsc2) / (I * omega - (e1 - esc2)) +
                               (wsc1 + w1) / (I * omega - (esc1 - e1)) -(wsc2 + w1) / (I * omega - (esc2 - e1)));
    };

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + Pairing("A", Delta);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);
    INFO("Energy levels " << H.getEigenValues());
    INFO("The value of ground energy is " << H.getGroundEnergy());

    DensityMatrix rho(S, H, beta);
    rho.prepare();
    rho.compute();
    for(QuantumState i = 0; i < S.getNumberOfStates(); ++i)
        INFO("Weight " << i << " = " << rho.getWeight(i));

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex up_index = IndexInfo.getIndex("A", 0, up);
    ParticleIndex down_index = IndexInfo.getIndex("A", 0, down);

    GreensFunction G(S,
                     H,
                     Operators.getAnnihilationOperator(down_index),
                     Operators.getCreationOperator(down_index),
                     rho);
    GreensFunction F(S,
                     H,
                     Operators.getAnnihilationOperator(up_index),
                     Operators.getAnnihilationOperator(down_index),
                     rho);

    G.prepare();
    G.compute();
    F.prepare();
    F.compute();

    for(int n = 0; n < 100; ++n) {
        REQUIRE_THAT(G(n), IsCloseTo(G_ref(n), 1e-14));
        REQUIRE_THAT(F(n), IsCloseTo(F_ref(n), 1e-14));
    }

    // cppcheck-suppress syntaxError
    SECTION("GFContainer") {
        GFContainer G(IndexInfo, S, H, rho, Operators);
        GFContainer F(IndexInfo, S, H, rho, Operators, true);

        std::set<IndexCombination2> indices;
        indices.insert(IndexCombination2(0, 0));
        indices.insert(IndexCombination2(0, 1));
        indices.insert(IndexCombination2(1, 0));
        indices.insert(IndexCombination2(1, 1));

        G.prepareAll(indices);
        G.computeAll();
        F.prepareAll(indices);
        F.computeAll();

        for(int n = -100; n < 100; ++n) {
            REQUIRE_THAT(G(0, 0)(n), IsCloseTo(G_ref(n), 1e-14));
            REQUIRE_THAT(G(0, 1)(n), IsCloseTo(0, 1e-14));
            REQUIRE_THAT(G(1, 0)(n), IsCloseTo(0, 1e-14));
            REQUIRE_THAT(G(1, 1)(n), IsCloseTo(G_ref(n), 1e-14));

            REQUIRE_THAT(F(0, 0)(n), IsCloseTo(0, 1e-14));
            REQUIRE_THAT(F(0, 1)(n), IsCloseTo(-F_ref(n), 1e-14));
            REQUIRE_THAT(F(1, 0)(n), IsCloseTo(F_ref(n), 1e-14));
            REQUIRE_THAT(F(1, 1)(n), IsCloseTo(0, 1e-14));
        }
    }
}
