//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <pomerol/DensityMatrix.hpp>
#include <pomerol/EnsembleAverage.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/MonomialOperator.hpp>
#include <pomerol/Susceptibility.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>

using namespace Pomerol;

TEST_CASE("Susceptibilities of a single Hubbard atom", "[Susceptibility]") {
    RealType U = 1.0;
    RealType mu = 0.4;
    RealType h_field = 0.01;
    RealType beta = 10.0;
    int n_iw = 20;

    auto omega = [beta](int n) { return M_PI * (2 * n) / beta; };

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + Magnetization("A", -h_field);
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

    ParticleIndex up_index = IndexInfo.getIndex("A", 0, up);
    ParticleIndex dn_index = IndexInfo.getIndex("A", 0, down);

    // Quadratic operators of form c^+ c
    QuadraticOperator s_plus(IndexInfo, HS, S, H, up_index, dn_index);
    QuadraticOperator s_minus(IndexInfo, HS, S, H, dn_index, up_index);
    QuadraticOperator n_up(IndexInfo, HS, S, H, up_index, up_index);
    QuadraticOperator n_dn(IndexInfo, HS, S, H, dn_index, dn_index);

    for(auto* op : {&s_plus, &s_minus, &n_up, &n_dn}) {
        op->prepare(HS);
        op->compute();
    }

    // Reference statistical weights of states
    RealVectorType weights_ref(4); // {0, up, down, 2}
    weights_ref << 1.0, exp(-beta * (-mu - h_field)), exp(-beta * (-mu + h_field)), exp(-beta * (-2 * mu + U));
    weights_ref /= weights_ref.sum();
    RealType wu = weights_ref[1];
    RealType wd = weights_ref[2];
    RealType w2 = weights_ref[3];

    SECTION("Ensemble averages") {
        EnsembleAverage s_plus_aver(S, H, s_plus, rho);
        s_plus_aver.compute();
        REQUIRE_THAT(s_plus_aver.getResult(), IsCloseTo(0, 1e-14));

        EnsembleAverage s_minus_aver(S, H, s_minus, rho);
        s_minus_aver.compute();
        REQUIRE_THAT(s_minus_aver.getResult(), IsCloseTo(0, 1e-14));

        EnsembleAverage n_up_aver(S, H, n_up, rho);
        n_up_aver.compute();
        RealType n_up_ref = weights_ref[1] + weights_ref[3];
        REQUIRE_THAT(n_up_aver.getResult(), IsCloseTo(n_up_ref, 1e-14));

        EnsembleAverage n_dn_aver(S, H, n_dn, rho);
        n_dn_aver.compute();
        RealType n_dn_ref = weights_ref[2] + weights_ref[3];
        REQUIRE_THAT(n_dn_aver.getResult(), IsCloseTo(n_dn_ref, 1e-14));
    }

    SECTION("<S_+; S_->") {
        Susceptibility Chi(S, H, s_plus, s_minus, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) {
            ComplexType g = 0;
            if(std::abs(wu - wd) < 1e-8 && n == 0) // E_up == E_down
                g += wu * beta;
            else
                g += -(wu - wd) / (I * omega(n) - 2 * h_field);
            return g;
        };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }

    SECTION("<n_up; n_up>") {
        Susceptibility Chi(S, H, n_up, n_up, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) { return n == 0 ? (wu + w2) * (1 - wu - w2) * beta : 0; };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }

    SECTION("<n_up; n_dn>") {
        Susceptibility Chi(S, H, n_up, n_dn, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) { return n == 0 ? (w2 - (wu + w2) * (wd + w2)) * beta : 0; };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }
}
