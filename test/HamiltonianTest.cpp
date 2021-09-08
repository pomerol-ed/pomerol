//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2021 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** \file tests/IndexTerm.cpp
** \brief Test of the Symmetrizer::IndexPermutation.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include <pomerol/Hamiltonian.hpp>
#include <pomerol/HamiltonianPart.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/Index.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/Misc.hpp>
#include <pomerol/MonomialOperator.hpp>
#include <pomerol/MonomialOperatorPart.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/StatesClassification.hpp>

#include <mpi_dispatcher/misc.hpp>

#include "catch2/catch-pomerol.hpp"

using namespace Pomerol;

TEST_CASE("Simple Hamiltonian test", "[hamiltonian]") {
    using namespace LatticePresets;

    auto HExpr = CoulombS("A", 1.0, -0.5);
    HExpr += CoulombS("B", 2.0, -1.0);
    HExpr += Hopping("A", "B", -1.0);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO(IndexInfo);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    INFO(H.getPart(BlockNumber(4)));
    H.compute(MPI_COMM_WORLD);
    INFO(H.getPart(BlockNumber(4)));

    SECTION("Ground state energy") {
        RealType E_ref = -2.8860009;
        RealType E = H.getGroundEnergy();

        REQUIRE_THAT(E, IsCloseTo(E_ref, 1e-7));
    }

    SECTION("Monomial operators") {
        ParticleIndex op_index = IndexInfo.getIndex("B", 0, up);
        BlockNumber test_block = 4;

        CreationOperator op(IndexInfo, HS, S, H, op_index);
        op.prepare(HS);
        op.compute();
        BlockNumber result_block = op.getLeftIndex(test_block);

        INFO("Acting with rotated cdag_" << op_index << " on block " << test_block << " and receiving "
                                         << result_block);

        using LOperatorT = LOperatorType<RealType>;

        HamiltonianPart HpartRHS(LOperatorT(HExpr, HS.getFullHilbertSpace()), S, test_block);
        HpartRHS.prepare();
        HpartRHS.compute();
        INFO(HpartRHS);

        HamiltonianPart HpartLHS(LOperatorT(HExpr, HS.getFullHilbertSpace()), S, result_block);
        HpartLHS.prepare();
        HpartLHS.compute();
        INFO(HpartLHS);

        auto cdag1op = LOperatorT(Operators::c_dag("B", (unsigned short)0, up), HS.getFullHilbertSpace());
        MonomialOperatorPart Cdag1(cdag1op, S, HpartRHS, HpartLHS);
        Cdag1.compute();
        INFO(Cdag1);

        auto c1op = LOperatorT(Operators::c("B", (unsigned short)0, up), HS.getFullHilbertSpace());
        MonomialOperatorPart C1(c1op, S, HpartLHS, HpartRHS);
        C1.compute();
        INFO(C1);

        // Check transposition
        auto diff1 = (Cdag1.getRowMajorValue<false>() - C1.getColMajorValue<false>().transpose()).eval();
        diff1.prune(1e-12);
        REQUIRE(diff1.nonZeros() == 0);

        auto diff2 = (Cdag1.getColMajorValue<false>() - C1.getRowMajorValue<false>().transpose()).eval();
        diff2.prune(1e-12);
        REQUIRE(diff2.nonZeros() == 0);
    }
}
