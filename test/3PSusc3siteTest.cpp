//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/3PSusc3siteTest.cpp
/// \brief 3-point susceptibilities of a small Hubbard cluster.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/MonomialOperator.hpp>
#include <pomerol/ThreePointSusceptibility.hpp>

#include "catch2/catch-pomerol.hpp"

using namespace Pomerol;

TEST_CASE("3-point susceptibilities of a small Hubbard cluster", "[ThreePointSusceptibility]") {
    RealType U = 4.0;
    RealType mu = 0.6 * U;
    RealType t = 1.0;
    RealType beta = 5.0;

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + CoulombS("B", U, -mu) + CoulombS("C", U, -mu);
    HExpr += Hopping("A", "B", -t);
    HExpr += Hopping("B", "C", -t);
    HExpr += Hopping("C", "A", -t);
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

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex A_up_index = IndexInfo.getIndex("A", 0, up);
    ParticleIndex A_dn_index = IndexInfo.getIndex("A", 0, down);
    ParticleIndex C_up_index = IndexInfo.getIndex("C", 0, up);
    ParticleIndex C_dn_index = IndexInfo.getIndex("C", 0, down);

    SECTION("Particle-particle channel") {
        // Pairing operator
        QuadraticOperator Delta(IndexInfo, HS, S, H, C_up_index, C_dn_index, {false, false});
        Delta.prepare(HS);
        Delta.compute();

        ThreePointSusceptibility chi3pp(S,
                                        H,
                                        Operators.getCreationOperator(A_up_index),
                                        Operators.getCreationOperator(A_dn_index),
                                        Delta,
                                        rho);
        chi3pp.prepare();
        chi3pp.compute();

        // Reference values from 'chi3cluster.py'
        ComplexMatrixType chi3_ref(3, 3);
        chi3_ref << 0.056123393380680404 - 0.027886788783494966 * I, 0.03544036019820439,
            0.008008339386258077 - 0.0022610905291281 * I, 0.03544036019820439,
            0.056123393380680404 + 0.027886788783494966 * I, 0.026205160724877802 + 0.020851328654297216 * I,
            0.008008339386258187 - 0.002261090529127845 * I, 0.026205160724877885 + 0.020851328654297483 * I,
            0.013637220898198006 + 0.012985397056265282 * I;

        for(int n1 = -1; n1 <= 1; ++n1) {
            for(int n2 = -1; n2 <= 1; ++n2) {
                auto result = chi3pp(n1, n2);
                auto ref = chi3_ref(n1 + 1, n2 + 1);
                REQUIRE_THAT(result, IsCloseTo(ref, 1e-10));
            }
        }
    }

    SECTION("Particle-hole channel") {
        // Density operator
        QuadraticOperator N(IndexInfo, HS, S, H, C_dn_index, C_dn_index);
        N.prepare(HS);
        N.compute();

        ThreePointSusceptibility chi3ph(S,
                                        H,
                                        Operators.getCreationOperator(A_up_index),
                                        Operators.getAnnihilationOperator(A_up_index),
                                        N,
                                        rho);
        chi3ph.prepare();
        chi3ph.compute();

        // Reference values from 'chi3cluster.py'
        ComplexMatrixType chi3_ref(3, 3);
        chi3_ref << -0.03377687240880357 + 0.6799848693219751 * I, -0.01589165211289571,
            -0.010498942766103334 - 0.0019558762969712905 * I, -0.01589165211289571,
            -0.03377687240880357 - 0.6799848693219751 * I, 0.028277168735633868 + 0.006544223759963227 * I,
            -0.010498942766103213 - 0.001955876296971393 * I, 0.02827716873563382 + 0.006544223759962816 * I,
            0.03857021454863498 - 0.6687150248684948 * I;

        for(int n1 = -1; n1 <= 1; ++n1) {
            for(int n2 = -1; n2 <= 1; ++n2) {
                auto result = chi3ph(n1, n2);
                auto ref = chi3_ref(n1 + 1, n2 + 1);
                REQUIRE_THAT(result, IsCloseTo(ref, 1e-10));
            }
        }
    }
}
