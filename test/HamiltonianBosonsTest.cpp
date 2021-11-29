//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/HamiltonianBosonsTest.cpp
/// \brief Diagonalization of a Hubbard atom coupled to a boson.
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/StatesClassification.hpp>

#include <mpi_dispatcher/misc.hpp>

#include "catch2/catch-pomerol.hpp"

#include <algorithm>
#include <map>
#include <vector>

using namespace Pomerol;

TEST_CASE("Hamiltonian of an isolated Hubbard-Holstein atom", "[HubbardHolstein]") {
    using namespace LatticePresets;
    using namespace Operators;

    RealType U = 1.0;
    RealType mu = 0.4;
    RealType Omega = 4.0;
    RealType lambda = 1.5;
    unsigned int bits_per_boson = 6;
    unsigned int n_ev_to_check = 64;

    //
    // Reference eigenvalues
    //
    // NB.: These are exact eigenvalues of the full infinitely-dimensional
    // problem. As we truncate the bosonic Hilbert space, energies of the highly
    // excited states will deviate from the exact values. That's why we check
    // only the lowest n_ev_to_check eigenvalues.
    //

    // Renormalized chemical potential and Coulomb interaction
    RealType mu_r = mu + lambda * lambda / Omega;
    RealType U_r = U - 2 * lambda * lambda / Omega;

    std::vector<double> ev_ref;
    for(int nb = 0; nb < (1 << bits_per_boson); ++nb) {
        ev_ref.push_back(Omega * nb);
        ev_ref.push_back(-mu_r + Omega * nb);
        ev_ref.push_back(-mu_r + Omega * nb);
        ev_ref.push_back(-2 * mu_r + U_r + Omega * nb);
    }
    std::sort(ev_ref.begin(), ev_ref.end());

    auto HExpr = CoulombS("A", U, -mu);
    HExpr += BosonLevel("A", Omega, 0) + HolsteinInteraction("A", lambda, 0, 0);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    SECTION("bits_per_boson") {
        auto HS = MakeHilbertSpace(IndexInfo, HExpr, bits_per_boson);
        HS.compute();
        StatesClassification S;
        S.compute(HS);

        Hamiltonian H(S);
        H.prepare(HExpr, HS, MPI_COMM_WORLD);
        H.compute(MPI_COMM_WORLD);

        // Sorted list of eigenvalues
        auto ev = H.getEigenValues();
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        std::sort(ev.data(), ev.data() + ev.size());
        INFO(ev);

        REQUIRE(ev.size() == ev_ref.size());
        for(int n = 0; n < n_ev_to_check; ++n)
            REQUIRE_THAT(ev(n), IsCloseTo(ev_ref[n], 1e-10));
    }

    SECTION("bits_per_boson_map") {
        std::map<decltype(HExpr)::index_types, unsigned int> bits_per_boson_map;
        auto boson_indices = decltype(HExpr)::index_types("A", 0, undef);
        bits_per_boson_map[boson_indices] = bits_per_boson;
        auto HS = MakeHilbertSpace(IndexInfo, HExpr, bits_per_boson_map);
        HS.compute();
        StatesClassification S;
        S.compute(HS);

        Hamiltonian H(S);
        H.prepare(HExpr, HS, MPI_COMM_WORLD);
        H.compute(MPI_COMM_WORLD);

        // Sorted list of eigenvalues
        auto ev = H.getEigenValues();
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        std::sort(ev.data(), ev.data() + ev.size());
        INFO(ev);

        REQUIRE(ev.size() == ev_ref.size());
        for(int n = 0; n < n_ev_to_check; ++n)
            REQUIRE_THAT(ev(n), IsCloseTo(ev_ref[n], 1e-10));
    }
}
