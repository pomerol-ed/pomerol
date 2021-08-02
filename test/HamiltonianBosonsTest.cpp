//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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

/** \file tests/green.cpp
** \brief Test of a Green's function calculation (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
*/

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

TEST_CASE("Hamiltonian of an isolated Hubbard-Holstein atom",
          "[HubbardHolstein]") {
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
    HExpr += BosonLevel("A", Omega, 0) +
             HolsteinInteraction("A", lambda, 0, 0);
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
        std::sort(ev.data(), ev.data() + ev.size());
        INFO(ev);

        REQUIRE(ev.size() == ev_ref.size());
        for(int n = 0; n < n_ev_to_check; ++n)
            REQUIRE_THAT(ev(n), IsCloseTo(ev_ref[n], 1e-10));
    }
}
