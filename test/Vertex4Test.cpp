//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/Vertex4Test.cpp
/// \brief Two-particle Green's function and the vertex of a single Hubbard atom.
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)
/// \author Igor Krivenko (igor.s.krivenko@gmail.com)

#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GFContainer.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/TwoParticleGFContainer.hpp>
#include <pomerol/Vertex4.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>

// Generalized 'square' function.
template <typename T> inline T sqr(T x) {
    return x * x;
}

using namespace Pomerol;

TEST_CASE("Two-particle vertex of a single Hubbard atom", "[Vertex4]") {
    // Parameters
    RealType U = 1.0;
    RealType beta = 40;
    int n_iw = 4;

    // Auxiliary functions
    auto w = [beta](int n) { return M_PI * (2 * n + 1) / beta; };
    auto delta = [](int n1, int n2) { return n1 == n2 ? 1.0 : 0.0; };
    auto deltam = [](int n1, int n2) { return n1 + n2 == -1 ? 1.0 : 0.0; };

    // Reference gamma4(up,up,up,up)
    auto Gamma4uuuu_ref = [&](int n1, int n2, int n3) {
        RealType omega1 = w(n1);
        RealType omega2 = w(n2);

        return -beta * (delta(n1, n3) - delta(n2, n3)) * sqr(0.5 * U) * (1. + sqr(0.5 * U / omega1)) *
               (1. + sqr(0.5 * U / omega2));
    };

    // Reference gamma4(up,down,up,down)
    auto Gamma4udud_ref = [&](int n1, int n2, int n3) {
        RealType omega1 = w(n1);
        RealType omega2 = w(n2);
        RealType omega3 = w(n3);
        RealType omega4 = omega1 + omega2 - omega3;

        RealType w = 1.0 / (1.0 + std::exp(beta * 0.5 * U));

        ComplexType Value = U;
        Value += -0.125 * U * U * U * (sqr(omega1) + sqr(omega2) + sqr(omega3) + sqr(omega4)) /
                 (omega1 * omega2 * omega3 * omega4);
        Value += -0.1875 * U * U * U * U * U / (omega1 * omega2 * omega3 * omega4);
        Value += -beta * (2 * deltam(n1, n2) + delta(n1, n3)) * w * sqr(0.5 * U) * (1.0 + sqr(0.5 * U / omega2)) *
                 (1.0 + sqr(0.5 * U / omega3));
        Value += beta * (2 * delta(n2, n3) + delta(n1, n3)) * (1 - w) * sqr(0.5 * U) * (1.0 + sqr(0.5 * U / omega1)) *
                 (1.0 + sqr(0.5 * U / omega2));

        return Value;
    };

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -U / 2);
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

    GFContainer G(IndexInfo, S, H, rho, Operators);
    G.prepareAll();
    G.computeAll();

    TwoParticleGFContainer Chi(IndexInfo, S, H, rho, Operators);
    Chi.ReduceResonanceTolerance = 1e-4;
    Chi.CoefficientTolerance = 1e-12;
    Chi.prepareAll();
    Chi.computeAll();
    MPI_Barrier(MPI_COMM_WORLD);

    ParticleIndex up_index = IndexInfo.getIndex("A", 0, up);
    ParticleIndex down_index = IndexInfo.getIndex("A", 0, down);

    SECTION("\\chi_{\\up\\up\\up\\up} and \\Gamma_{\\up\\up\\up\\up}") {
        GreensFunction const& GF = G(up_index, up_index);
        TwoParticleGF const& Chi_uuuu = Chi(IndexCombination4(up_index, up_index, up_index, up_index));
        Vertex4 Gamma4_uuuu(Chi_uuuu, GF, GF, GF, GF);
        Gamma4_uuuu.compute(n_iw);

        for(int n1 = -n_iw; n1 < n_iw; ++n1) {
            for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                for(int n3 = -n_iw; n3 < n_iw; ++n3) {
                    int n4 = n1 + n2 - n3;
                    INFO(n1 << " " << n2 << " " << n3 << " " << n1 + n2 - n3);

                    ComplexType chi_value = Chi_uuuu(n1, n2, n3);
                    ComplexType chi_ref = Gamma4uuuu_ref(n1, n2, n3) * GF(n1) * GF(n2) * GF(n3) * GF(n4) +
                                          beta * GF(n1) * GF(n2) * delta(n1, n4) -
                                          beta * GF(n1) * GF(n2) * delta(n1, n3);
                    REQUIRE_THAT(chi_value, IsCloseTo(chi_ref, 1e-6));

                    ComplexType Gamma_value = Gamma4_uuuu(n1, n2, n3);
                    ComplexType Gamma_ref = Gamma4uuuu_ref(n1, n2, n3) * GF(n1) * GF(n2) * GF(n3) * GF(n4);
                    REQUIRE_THAT(Gamma_value, IsCloseTo(Gamma_ref, 1e-10));
                }
            }
        }
    }

    SECTION("\\chi_{\\up\\down\\up\\down}") {
        GreensFunction const& GF_up = G(up_index, up_index);
        GreensFunction const& GF_down = G(down_index, down_index);
        TwoParticleGF const& Chi_udud = Chi(IndexCombination4(up_index, down_index, up_index, down_index));
        for(int n1 = -n_iw; n1 < n_iw; ++n1) {
            for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                for(int n3 = -n_iw; n3 < n_iw; ++n3) {
                    int n4 = n1 + n2 - n3;
                    INFO(n1 << " " << n2 << " " << n3 << " " << n4);

                    ComplexType chi_value = Chi_udud(n1, n2, n3);
                    ComplexType chi_ref =
                        Gamma4udud_ref(n1, n2, n3) * GF_up(n1) * GF_down(n2) * GF_up(n3) * GF_down(n4) -
                        beta * GF_up(n1) * GF_down(n2) * delta(n1, n3);
                    REQUIRE_THAT(chi_value, IsCloseTo(chi_ref, 1e-10));
                }
            }
        }
    }
}
