//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/GF2siteTest.cpp
/// \brief Single-particle Green's functions of a Hubbard dimer.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko

#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GreensFunction.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/StatesClassification.hpp>

#include "catch2/catch-pomerol.hpp"

#include <vector>

using namespace Pomerol;

TEST_CASE("Green's function of a Hubbard dimer", "[GF2site]") {
    RealType U = 1.0;
    RealType mu = 0.5;
    RealType beta = 10.0;

    // Reference Green's function
    ComplexVectorType G_ref(10);
    // clang-format off
    G_ref << -2.53021005e-01 * I,
             -4.62090702e-01 * I,
             -4.32482782e-01 * I,
             -3.65598615e-01 * I,
             -3.07785174e-01 * I,
             -2.62894141e-01 * I,
             -2.28274316e-01 * I,
             -2.01170772e-01 * I,
             -1.79539602e-01 * I, // cppcheck-suppress constStatement
             -1.61950993e-01 * I;
    // clang-format on

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + CoulombS("B", U, -mu);
    HExpr += Hopping("A", "B", -1.0);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    auto const& A_dn = IndexInfo.getInfo(0);
    auto const& A_up = IndexInfo.getInfo(1);
    auto const& B_dn = IndexInfo.getInfo(2);
    auto const& B_up = IndexInfo.getInfo(3);

    auto N = Operators::N(std::vector<decltype(IndexInfo)::IndexInfo>{A_up, A_dn, B_up, B_dn});
    INFO("N = " << N);

    auto Sz = Operators::Sz(std::vector<decltype(IndexInfo)::IndexInfo>{A_up, B_up},
                            std::vector<decltype(IndexInfo)::IndexInfo>{A_dn, B_dn});
    INFO("Sz = " << Sz);

    REQUIRE(Sz * N == N * Sz);
    REQUIRE(HExpr * N == N * HExpr);
    REQUIRE(HExpr * Sz == Sz * HExpr);

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

    ParticleIndex A_down_index = IndexInfo.getIndex("A", 0, down);

    auto const& c_map = Operators.getCreationOperator(A_down_index).getBlockMapping();
    for(auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++) {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(A_down_index),
                      Operators.getCreationOperator(A_down_index),
                      rho);

    GF.prepare();
    GF.compute();

    for(int n = 0; n < G_ref.size(); ++n) {
        auto result = GF(n);
        auto ref = G_ref[n];
        // The tolerance has to be fairly large as Pomerol discards
        // some contributions to the GF.
        REQUIRE_THAT(result, IsCloseTo(ref, 1e-8));
    }
}
