//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2025 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/ThreePointSusceptibilityPart.cpp
/// \brief Part of a 3-point susceptibility in the Matsubara representation (implementation).
/// \author Igor Krivenko

#include "pomerol/ThreePointSusceptibilityPart.hpp"
#include "pomerol/ChaseIndices.hpp"

#include <cassert>
#include <mutex>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex NonResonantFFTerm_mpi_datatype_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex NonResonantFBTerm_mpi_datatype_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex ResonantTerm_mpi_datatype_mutex;

namespace Pomerol {

//
// ThreePointSusceptibilityPart::NonResonantFFTerm
//

ThreePointSusceptibilityPart::NonResonantFFTerm&
ThreePointSusceptibilityPart::NonResonantFFTerm::operator+=(NonResonantFFTerm const& AnotherTerm) {
    long combinedWeight = Weight + AnotherTerm.Weight;
    for(unsigned short p = 0; p < 2; ++p)
        // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
        Poles[p] = (Weight * Poles[p] + AnotherTerm.Weight * AnotherTerm.Poles[p]) / combinedWeight;
    Weight = combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

// When called for the first time, this function creates an MPI structure datatype
// describing ThreePointSusceptibilityPart::NonResonantFFTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
MPI_Datatype ThreePointSusceptibilityPart::NonResonantFFTerm::mpi_datatype() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    std::lock_guard<std::mutex> const lock(NonResonantFFTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        int blocklengths[] = {1, 2, 1};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Aint displacements[] = {offsetof(NonResonantFFTerm, Coeff),
                                    offsetof(NonResonantFFTerm, Poles),
                                    offsetof(NonResonantFFTerm, Weight)};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Datatype types[] = {
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType Coeff
            MPI_DOUBLE,                 // RealType Poles[2]
            MPI_LONG                    // long Weight
        };
        MPI_Type_create_struct(3, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

//
// ThreePointSusceptibilityPart::NonResonantFBTerm
//

ThreePointSusceptibilityPart::NonResonantFBTerm&
ThreePointSusceptibilityPart::NonResonantFBTerm::operator+=(NonResonantFBTerm const& AnotherTerm) {
    long combinedWeight = Weight + AnotherTerm.Weight;
    // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
    P1 = (Weight * P1 + AnotherTerm.Weight * AnotherTerm.P1) / combinedWeight;
    // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
    P12 = (Weight * P12 + AnotherTerm.Weight * AnotherTerm.P12) / combinedWeight;
    Weight = combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

// When called for the first time, this function creates an MPI structure datatype
// describing ThreePointSusceptibilityPart::NonResonantFBTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
MPI_Datatype ThreePointSusceptibilityPart::NonResonantFBTerm::mpi_datatype() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    std::lock_guard<std::mutex> const lock(NonResonantFBTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        int blocklengths[] = {1, 1, 1, 1, 1};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Aint displacements[] = {offsetof(NonResonantFBTerm, Coeff),
                                    offsetof(NonResonantFBTerm, P1),
                                    offsetof(NonResonantFBTerm, P12),
                                    offsetof(NonResonantFBTerm, xi),
                                    offsetof(NonResonantFBTerm, Weight)};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Datatype types[] = {
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType Coeff
            MPI_DOUBLE,                 // RealType P1
            MPI_DOUBLE,                 // RealType P12
            MPI_INT,                    // int xi
            MPI_LONG                    // long Weight
        };
        MPI_Type_create_struct(5, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

//
// ThreePointSusceptibilityPart::ResonantTerm
//

ThreePointSusceptibilityPart::ResonantTerm&
ThreePointSusceptibilityPart::ResonantTerm::operator+=(ResonantTerm const& AnotherTerm) {
    long combinedWeight = Weight + AnotherTerm.Weight;
    // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
    P = (Weight * P + AnotherTerm.Weight * AnotherTerm.P) / combinedWeight;
    Weight = combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

// When called for the first time, this function creates an MPI structure datatype
// describing ThreePointSusceptibilityPart::ResonantTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
MPI_Datatype ThreePointSusceptibilityPart::ResonantTerm::mpi_datatype() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    std::lock_guard<std::mutex> const lock(ResonantTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        int blocklengths[] = {1, 1, 1, 1};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Aint displacements[] = {offsetof(ResonantTerm, Coeff),
                                    offsetof(ResonantTerm, P),
                                    offsetof(ResonantTerm, xi),
                                    offsetof(ResonantTerm, Weight)};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Datatype types[] = {
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType Coeff
            MPI_DOUBLE,                 // RealType P
            MPI_INT,                    // int xi
            MPI_LONG                    // long Weight
        };
        MPI_Type_create_struct(4, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

//
// ThreePointSusceptibilityPart
//

ThreePointSusceptibilityPart::ThreePointSusceptibilityPart(MonomialOperatorPart const& F1,
                                                           MonomialOperatorPart const& F2,
                                                           MonomialOperatorPart const& B1,
                                                           MonomialOperatorPart const& B2,
                                                           HamiltonianPart const& Hpart1,
                                                           HamiltonianPart const& Hpart2,
                                                           HamiltonianPart const& Hpart3,
                                                           DensityMatrixPart const& DMpart1,
                                                           DensityMatrixPart const& DMpart2,
                                                           DensityMatrixPart const& DMpart3,
                                                           Channel channel,
                                                           bool SwappedFermionOps)
    : Thermal(DMpart1.beta),
      ComputableObject(),
      F1(F1),
      F2(F2),
      B1(B1),
      B2(B2),
      Hpart1(Hpart1),
      Hpart2(Hpart2),
      Hpart3(Hpart3),
      DMpart1(DMpart1),
      DMpart2(DMpart2),
      DMpart3(DMpart3),
      channel(channel),
      SwappedFermionOps(SwappedFermionOps),
      NonResonantFFTerms(NonResonantFFTerm::Hash(), NonResonantFFTerm::KeyEqual(), NonResonantFFTerm::IsNegligible()),
      NonResonantFBTerms(NonResonantFBTerm::Hash(), NonResonantFBTerm::KeyEqual(), NonResonantFBTerm::IsNegligible()),
      ResonantTerms(ResonantTerm::Hash(), ResonantTerm::KeyEqual(), ResonantTerm::IsNegligible()) {}

void ThreePointSusceptibilityPart::compute() {
    if(getStatus() >= Computed)
        return;

    if(F1.isComplex() || F2.isComplex() || B1.isComplex() || B2.isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();
}

template <bool Complex> void ThreePointSusceptibilityPart::computeImpl() {
    NonResonantFFTerms.clear();
    NonResonantFBTerms.clear();
    ResonantTerms.clear();

    RealType beta = DMpart1.beta;
    RealType prefactor = ((channel == PH) == SwappedFermionOps) ? 1 : -1;

    // B = B_1 B_2
    ColMajorMatrixType<Complex> Bmatrix = (B1.getColMajorValue<Complex>() * B2.getColMajorValue<Complex>()).pruned();

    // <1| F1 |2> <2| F2 |3> <3| B |1>
    RowMajorMatrixType<Complex> const& F1matrix = F1.getRowMajorValue<Complex>();
    RowMajorMatrixType<Complex> const& F2matrix = F2.getRowMajorValue<Complex>();

    InnerQuantumState index1Max = Bmatrix.outerSize();
    for(InnerQuantumState index1 = 0; index1 < index1Max; ++index1) {
        RealType E1 = Hpart1.getEigenValue(index1);
        RealType weight1 = DMpart1.getWeight(index1);

        typename RowMajorMatrixType<Complex>::InnerIterator index2ket_iter(F1matrix, index1);
        for(; index2ket_iter; ++index2ket_iter) {
            RealType E2 = Hpart2.getEigenValue(index2ket_iter.index());
            RealType weight2 = DMpart2.getWeight(index2ket_iter.index());

            typename ColMajorMatrixType<Complex>::InnerIterator index3bra_iter(Bmatrix, index1);
            typename RowMajorMatrixType<Complex>::InnerIterator index3ket_iter(F2matrix, index2ket_iter.index());
            while(index3bra_iter && index3ket_iter) {
                if(chaseIndices<Complex>(index3ket_iter, index3bra_iter)) {
                    RealType E3 = Hpart3.getEigenValue(index3ket_iter.index());
                    RealType weight3 = DMpart3.getWeight(index3ket_iter.index());

                    ComplexType MatrixElement =
                        index2ket_iter.value() * index3ket_iter.value() * index3bra_iter.value();
                    MatrixElement *= prefactor;

                    addMultiterm(MatrixElement, beta, E1, E2, E3, weight1, weight2, weight3);

                    ++index3bra_iter;
                    ++index3ket_iter;
                }
            }
        }
    }

    INFO("Total " << NonResonantFFTerms.size() << "+" << NonResonantFBTerms.size() << "+" << ResonantTerms.size() << "="
                  << NonResonantFFTerms.size() + NonResonantFBTerms.size() + ResonantTerms.size() << " terms");

    assert(NonResonantFFTerms.check_terms());
    assert(NonResonantFBTerms.check_terms());
    assert(ResonantTerms.check_terms());

    setStatus(Computed);
}

inline void ThreePointSusceptibilityPart::addMultiterm(ComplexType Coeff,
                                                       RealType beta,
                                                       RealType Ei,
                                                       RealType Ej,
                                                       RealType Ek,
                                                       RealType Wi,
                                                       RealType Wj,
                                                       RealType Wk) {
    RealType Eij = Ei - Ej;
    RealType Ejk = Ej - Ek;
    RealType Eik = Ei - Ek;

    // Non-resonant fermion-fermion term that is always present.
    ComplexType CoeffNRFF = Coeff * (Wi + Wj);
    if(std::abs(CoeffNRFF) > CoefficientTolerance) {
        NonResonantFFTerms.add_term(NonResonantFFTerm((channel == PP) ? -CoeffNRFF : CoeffNRFF,
                                                      SwappedFermionOps ? Ejk : Eij,
                                                      (SwappedFermionOps ? Eij : Ejk) * ((channel == PP) ? 1 : -1)));
    }

    // Resonant term
    if(std::abs(Eik) < ReduceResonanceTolerance) {
        ComplexType CoeffR = -Coeff * beta * Wi;
        if(std::abs(CoeffR) > CoefficientTolerance) {
            ResonantTerms.add_term(ResonantTerm(SwappedFermionOps ? -CoeffR : CoeffR,
                                                SwappedFermionOps ? Ejk : -Ejk,
                                                (channel == PP) ? -1 : 1));
        }
    } else {
        // Non-resonant fermion-boson term
        ComplexType CoeffNRFB = Coeff * (Wk - Wi);
        if(std::abs(CoeffNRFB) > CoefficientTolerance) {
            NonResonantFBTerms.add_term(NonResonantFBTerm(SwappedFermionOps ? -CoeffNRFB : CoeffNRFB,
                                                          SwappedFermionOps ? Ejk : Eij,
                                                          Eik,
                                                          (channel == PP) ? -1 : 1));

            // Extra non-resonant fermion-fermion term
            if(!SwappedFermionOps) {
                NonResonantFFTerms.add_term(
                    NonResonantFFTerm((channel == PP) ? -CoeffNRFB : CoeffNRFB, Eij, (channel == PP) ? Ejk : -Ejk));
            }
        }
    }
}

ComplexType ThreePointSusceptibilityPart::operator()(long MatsubaraNumber1, long MatsubaraNumber2) const {
    long MatsubaraNumberOdd1 = 2 * MatsubaraNumber1 + 1;
    long MatsubaraNumberOdd2 = 2 * MatsubaraNumber2 + 1;
    return (*this)(MatsubaraSpacing * RealType(MatsubaraNumberOdd1), MatsubaraSpacing * RealType(MatsubaraNumberOdd2));
}

ComplexType ThreePointSusceptibilityPart::operator()(ComplexType z1, ComplexType z2) const {
    if(getStatus() != Computed) {
        throw StatusMismatch("3PSusceptibilityPart: Calling operator() on uncomputed container. Did you purge all the "
                             "terms when called compute()?");
    }

    return NonResonantFFTerms(z1, z2) + NonResonantFBTerms(z1, z2) + ResonantTerms(z1, z2, ReduceResonanceTolerance);
}

void ThreePointSusceptibilityPart::clear() {
    NonResonantFFTerms.clear();
    NonResonantFBTerms.clear();
    ResonantTerms.clear();
    setStatus(Constructed);
}

} // namespace Pomerol
