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

GreensFunctionPart::GreensFunctionPart(MonomialOperatorPart const& C,
                                       MonomialOperatorPart const& CX,
                                       HamiltonianPart const& HpartInner,
                                       HamiltonianPart const& HpartOuter,
                                       DensityMatrixPart const& DMpartInner,
                                       DensityMatrixPart const& DMpartOuter)
    : Thermal(DMpartInner.beta),
      HpartInner(HpartInner),
      HpartOuter(HpartOuter),
      DMpartInner(DMpartInner),
      DMpartOuter(DMpartOuter),
      C(C),
      CX(CX),
      Terms(Term::Hash(), Term::KeyEqual(), Term::IsNegligible()) {}

void GreensFunctionPart::compute() {
    if(C.isComplex() || CX.isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();
}

template <bool Complex> void GreensFunctionPart::computeImpl() {
    Terms.clear();

    // Blocks (submatrices) of C and CX
    RowMajorMatrixType<Complex> const& Cmatrix = C.template getRowMajorValue<Complex>();
    ColMajorMatrixType<Complex> const& CXmatrix = CX.template getColMajorValue<Complex>();
    QuantumState outerSize = Cmatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Cmatrix and empty columns of CXmatrix.
    for(QuantumState index1 = 0; index1 < outerSize; ++index1) {
        // <index1|C|Cinner><CXinner|CX|index1>
        typename RowMajorMatrixType<Complex>::InnerIterator Cinner(Cmatrix, index1);
        typename ColMajorMatrixType<Complex>::InnerIterator CXinner(CXmatrix, index1);

        // While we are not at the last column of Cmatrix or at the last row of CXmatrix.
        while(Cinner && CXinner) {
            QuantumState C_index2 = Cinner.index();
            QuantumState CX_index2 = CXinner.index();

            // A meaningful matrix element
            if(C_index2 == CX_index2) {
                ComplexType Residue = Cinner.value() * CXinner.value() *
                                      (DMpartOuter.getWeight(index1) + DMpartInner.getWeight(C_index2));
                if(std::abs(Residue) > CoefficientTolerance) // Is the residue relevant?
                {
                    // Create a new term and append it to the list.
                    RealType Pole = HpartInner.getEigenValue(C_index2) - HpartOuter.getEigenValue(index1);
                    Terms.add_term(Term(Residue, Pole));
                };
                ++Cinner;  // The next non-zero element
                ++CXinner; // The next non-zero element
            } else {
                // Chasing: one index runs down the other index
                if(CX_index2 < C_index2)
                    for(; QuantumState(CXinner.index()) < C_index2; ++CXinner)
                        ;
                else
                    for(; QuantumState(Cinner.index()) < CX_index2; ++Cinner)
                        ;
            }
        }
    }

    assert(Terms.check_terms());
}

} // namespace Pomerol
