//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/GreensFunctionPart.cpp
/// \brief Part of a single-particle Matsubara Green's function (implementation).
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/GreensFunctionPart.hpp"

#include <cassert>
#include <cmath>

namespace Pomerol {

GreensFunctionPart::Term::Term(ComplexType Residue, RealType Pole) : Residue(Residue), Pole(Pole){};
ComplexType GreensFunctionPart::Term::operator()(ComplexType Frequency) const {
    return Residue / (Frequency - Pole);
}

ComplexType GreensFunctionPart::Term::operator()(RealType tau, RealType beta) const {
    using std::exp;
    return Pole > 0 ? -Residue * exp(-tau * Pole) / (1 + exp(-beta * Pole)) :
                      -Residue * exp((beta - tau) * Pole) / (exp(beta * Pole) + 1);
}

inline GreensFunctionPart::Term& GreensFunctionPart::Term::operator+=(Term const& AnotherTerm) {
    Residue += AnotherTerm.Residue;
    return *this;
}

GreensFunctionPart::GreensFunctionPart(MonomialOperatorPart const& F1,
                                       MonomialOperatorPart const& F2,
                                       HamiltonianPart const& HpartInner,
                                       HamiltonianPart const& HpartOuter,
                                       DensityMatrixPart const& DMpartInner,
                                       DensityMatrixPart const& DMpartOuter,
                                       RealType PoleResolution,
                                       RealType CoefficientTolerance)
    : Thermal(DMpartInner.beta),
      HpartInner(HpartInner),
      HpartOuter(HpartOuter),
      DMpartInner(DMpartInner),
      DMpartOuter(DMpartOuter),
      F1(F1),
      F2(F2),
      Terms(Term::Hash(PoleResolution), Term::KeyEqual(PoleResolution), Term::IsNegligible(CoefficientTolerance)),
      PoleResolution(PoleResolution),
      CoefficientTolerance(CoefficientTolerance) {}

void GreensFunctionPart::compute() {
    if(F1.isComplex() || F2.isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();
}

template <bool Complex> void GreensFunctionPart::computeImpl() {
    Terms.clear();

    // Blocks (submatrices) of F1 and F2
    RowMajorMatrixType<Complex> const& F1matrix = F1.template getRowMajorValue<Complex>();
    ColMajorMatrixType<Complex> const& F2matrix = F2.template getColMajorValue<Complex>();
    QuantumState outerSize = F1matrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of F1matrix and empty columns of F2matrix.
    for(QuantumState index1 = 0; index1 < outerSize; ++index1) {
        // <index1|F1|F1inner><F2inner|F2|index1>
        typename RowMajorMatrixType<Complex>::InnerIterator F1inner(F1matrix, index1);
        typename ColMajorMatrixType<Complex>::InnerIterator F2inner(F2matrix, index1);

        // While we are not at the last column of F1matrix or at the last row of F2matrix.
        while(F1inner && F2inner) {
            QuantumState F1_index2 = F1inner.index();
            QuantumState F2_index2 = F2inner.index();

            // A meaningful matrix element
            if(F1_index2 == F2_index2) {
                ComplexType Residue = F1inner.value() * F2inner.value() *
                                      (DMpartOuter.getWeight(index1) + DMpartInner.getWeight(F1_index2));
                if(std::abs(Residue) > CoefficientTolerance) // Is the residue relevant?
                {
                    // Create a new term and append it to the list.
                    RealType Pole = HpartInner.getEigenValue(F1_index2) - HpartOuter.getEigenValue(index1);
                    Terms.add_term(Term(Residue, Pole));
                };
                ++F1inner; // The next non-zero element
                ++F2inner; // The next non-zero element
            } else {
                // Chasing: one index runs down the other index
                if(F2_index2 < F1_index2)
                    for(; QuantumState(F2inner.index()) < F1_index2; ++F2inner)
                        ;
                else
                    for(; QuantumState(F1inner.index()) < F2_index2; ++F1inner)
                        ;
            }
        }
    }

    assert(Terms.check_terms());
}

} // namespace Pomerol
