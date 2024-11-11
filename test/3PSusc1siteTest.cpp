//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file test/3PSusc1siteTest.cpp
/// \brief 3-point susceptibilities of a single Hubbard atom.
/// \author Igor Krivenko

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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Pomerol;

// Auxiliary function g used in expressions for reference solutions
struct g_aux {
    RealType beta;
    // Order of many-body states: |0>, |up>, |dn>, |2>
    std::vector<RealType> E;
    std::vector<RealType> rho;

    g_aux(RealType beta, RealType U, RealType mu, RealType h_field)
        : beta(beta), E{0, -mu - h_field, -mu + h_field, U - 2 * mu}, rho(4) {
        std::transform(E.begin(), E.end(), rho.begin(), [beta](RealType e) { return std::exp(-beta * e); });
        double Z = std::accumulate(rho.begin(), rho.end(), .0);
        std::transform(rho.begin(), rho.end(), rho.begin(), [Z](RealType w) { return w / Z; });
    }

    static double delta(RealType w1, RealType w2) { return std::abs(w1 - w2) < 1e-14 ? 1.0 : 0.0; };

    ComplexType operator()(int i, int j, int k, RealType w1, RealType w2) const {
        if(std::abs(rho[k] - rho[i]) < 1e-14) {
            return 1. / (I * w2 + E[j] - E[k]) * (rho[i] + rho[j]) / (I * w1 + E[i] - E[j]) +
                   beta * delta(w1, -w2) * rho[i] / (I * w2 + E[j] - E[i]);
        } else {
            return 1. / (I * w2 + E[j] - E[k]) *
                   ((rho[k] - rho[i]) / (I * w1 + I * w2 + E[i] - E[k]) + (rho[i] + rho[j]) / (I * w1 + E[i] - E[j]));
        }
    }
};

TEST_CASE("3-point susceptibilities of a single Hubbard atom", "[ThreePointSusceptibility]") {
    RealType U = 1.0;
    RealType mu = 0.4;
    auto h_field = GENERATE(values({.0, 0.01}));
    RealType beta = 10.0;
    int n_iw = 20;

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

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex up_index = IndexInfo.getIndex("A", 0, up);
    ParticleIndex dn_index = IndexInfo.getIndex("A", 0, down);

    auto omega = [beta](int n) { return M_PI * (2 * n + 1) / beta; };
    g_aux g(beta, U, mu, h_field);

    // cppcheck-suppress syntaxError
    SECTION("Particle-particle channel") {
        for(auto index1 : {up_index, dn_index}) {
            for(auto index2 : {up_index, dn_index}) {
                ThreePointSusceptibility chi3pp(Channel::PP,
                                                S,
                                                H,
                                                Operators.getCreationOperator(index1),
                                                Operators.getAnnihilationOperator(index1),
                                                Operators.getCreationOperator(index2),
                                                Operators.getAnnihilationOperator(index2),
                                                rho);
                chi3pp.prepare();
                chi3pp.compute();

                if(index1 == index2) {
                    REQUIRE(chi3pp.isVanishing());
                } else {
                    auto ref = [&](int n1, int n2) {
                        return g(3, (index2 == up_index ? 1 : 2), 0, -omega(n1), -omega(n2)) +
                               g(3, (index1 == up_index ? 1 : 2), 0, -omega(n2), -omega(n1));
                    };
                    for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                        for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                            REQUIRE_THAT(chi3pp(n1, n2), IsCloseTo(ref(n1, n2), 1e-14));
                        }
                    }
                }
            }
        }
    }

    SECTION("Particle-hole channel") {
        for(auto index1 : {up_index, dn_index}) {
            for(auto index2 : {up_index, dn_index}) {
                ThreePointSusceptibility chi3ph(Channel::PH,
                                                S,
                                                H,
                                                Operators.getCreationOperator(index1),
                                                Operators.getAnnihilationOperator(index1),
                                                Operators.getCreationOperator(index2),
                                                Operators.getAnnihilationOperator(index2),
                                                rho);
                chi3ph.prepare();
                chi3ph.compute();

                if(index1 == index2) {
                    auto ref = [&](int n1, int n2) {
                        int st1 = index1 == up_index ? 1 : 2;
                        int st2 = index1 == up_index ? 2 : 1;
                        return g(st1, 0, st1, -omega(n1), omega(n2)) + g(3, st2, 3, -omega(n1), omega(n2));
                    };
                    for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                        for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                            REQUIRE_THAT(chi3ph(n1, n2), IsCloseTo(ref(n1, n2), 1e-14));
                        }
                    }
                } else {
                    auto ref = [&](int n1, int n2) {
                        int st = index2 == up_index ? 1 : 2;
                        return g(3, st, 3, -omega(n1), omega(n2)) - g(st, 3, st, omega(n1), -omega(n2));
                    };
                    for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                        for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                            REQUIRE_THAT(chi3ph(n1, n2), IsCloseTo(ref(n1, n2), 1e-14));
                        }
                    }
                }
            }
        }
    }

    SECTION("Crossed particle-hole channel") {
        for(auto index1 : {up_index, dn_index}) {
            for(auto index2 : {up_index, dn_index}) {
                ThreePointSusceptibility chi3xph(Channel::xPH,
                                                 S,
                                                 H,
                                                 Operators.getCreationOperator(index1),
                                                 Operators.getAnnihilationOperator(index1),
                                                 Operators.getCreationOperator(index2),
                                                 Operators.getAnnihilationOperator(index2),
                                                 rho);
                chi3xph.prepare();
                chi3xph.compute();

                if(index1 == index2) {
                    auto ref = [&](int n1, int n2) {
                        int st1 = index1 == up_index ? 1 : 2;
                        int st2 = index1 == up_index ? 2 : 1;
                        return -g(st1, 0, st1, -omega(n1), omega(n2)) - g(3, st2, 3, -omega(n1), omega(n2));
                    };
                    for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                        for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                            REQUIRE_THAT(chi3xph(n1, n2), IsCloseTo(ref(n1, n2), 1e-14));
                        }
                    }
                } else {
                    auto ref = [&](int n1, int n2) {
                        int st1 = index1 == up_index ? 1 : 2;
                        int st2 = index1 == up_index ? 2 : 1;
                        return -g(st1, 0, st2, -omega(n1), omega(n2)) - g(st1, 3, st2, omega(n2), -omega(n1));
                    };
                    for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                        for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                            REQUIRE_THAT(chi3xph(n1, n2), IsCloseTo(ref(n1, n2), 1e-14));
                        }
                    }
                }
            }
        }
    }

    SECTION("Crossing symmetry") {
        for(auto index1 : {up_index, dn_index}) {
            for(auto index2 : {up_index, dn_index}) {
                for(auto index3 : {up_index, dn_index}) {
                    for(auto index4 : {up_index, dn_index}) {
                        ThreePointSusceptibility chi3ph(Channel::PH,
                                                        S,
                                                        H,
                                                        Operators.getCreationOperator(index1),
                                                        Operators.getAnnihilationOperator(index2),
                                                        Operators.getCreationOperator(index3),
                                                        Operators.getAnnihilationOperator(index4),
                                                        rho);
                        chi3ph.prepare();
                        chi3ph.compute();

                        ThreePointSusceptibility chi3xph(Channel::xPH,
                                                         S,
                                                         H,
                                                         Operators.getCreationOperator(index1),
                                                         Operators.getAnnihilationOperator(index4),
                                                         Operators.getCreationOperator(index3),
                                                         Operators.getAnnihilationOperator(index2),
                                                         rho);
                        chi3xph.prepare();
                        chi3xph.compute();

                        for(int n1 = -n_iw; n1 < n_iw; ++n1) {
                            for(int n2 = -n_iw; n2 < n_iw; ++n2) {
                                REQUIRE_THAT(chi3xph(n1, n2), IsCloseTo(-chi3ph(n1, n2), 1e-14));
                            }
                        }
                    }
                }
            }
        }
    }
}
