//
// This file is part of pomerol, an exact diagonalization library aimed at
// solving condensed matter models of interacting fermions.
//
// Copyright (C) 2016-2024 A. Antipov, I. Krivenko and contributors
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// \file src/pomerol/TwoParticleGFPart.cpp
/// \brief Part of a fermionic two-particle Matsubara Green's function (implementation).
/// \author Igor Krivenko
/// \author Andrey Antipov (andrey.e.antipov@gmail.com)

#include "pomerol/TwoParticleGFPart.hpp"
#include "pomerol/ChaseIndices.hpp"

#include <cassert>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex NonResonantTerm_mpi_datatype_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex ResonantTerm_mpi_datatype_mutex;

namespace Pomerol {

//
// TwoParticleGFPart::NonResonantTerm
//
TwoParticleGFPart::NonResonantTerm& TwoParticleGFPart::NonResonantTerm::operator+=(NonResonantTerm const& AnotherTerm) {
    long combinedWeight = Weight + AnotherTerm.Weight;
    for(unsigned short p = 0; p < 3; ++p)
        // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
        Poles[p] = (Weight * Poles[p] + AnotherTerm.Weight * AnotherTerm.Poles[p]) / combinedWeight;
    Weight = combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

// When called for the first time, this function creates an MPI structure datatype
// describing TwoParticleGFPart::NonResonantTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
MPI_Datatype TwoParticleGFPart::NonResonantTerm::mpi_datatype() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    std::lock_guard<std::mutex> const lock(NonResonantTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        int blocklengths[] = {1, 3, 1, 1};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Aint displacements[] = {offsetof(NonResonantTerm, Coeff),
                                    offsetof(NonResonantTerm, Poles),
                                    offsetof(NonResonantTerm, isz4),
                                    offsetof(NonResonantTerm, Weight)};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Datatype types[] = {
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType Coeff
            MPI_DOUBLE,                 // RealType Poles[3]
            POMEROL_MPI_BOOL,           // bool isz4
            MPI_LONG                    // long Weight
        };
        MPI_Type_create_struct(4, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

//
// TwoParticleGFPart::ResonantTerm
//
TwoParticleGFPart::ResonantTerm& TwoParticleGFPart::ResonantTerm::operator+=(ResonantTerm const& AnotherTerm) {
    long combinedWeight = Weight + AnotherTerm.Weight;
    for(unsigned short p = 0; p < 3; ++p)
        // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
        Poles[p] = (Weight * Poles[p] + AnotherTerm.Weight * AnotherTerm.Poles[p]) / combinedWeight;
    Weight = combinedWeight;
    ResCoeff += AnotherTerm.ResCoeff;
    NonResCoeff += AnotherTerm.NonResCoeff;
    return *this;
}

// When called for the first time, this function creates an MPI structure datatype
// describing TwoParticleGFPart::ResonantTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
MPI_Datatype TwoParticleGFPart::ResonantTerm::mpi_datatype() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    std::lock_guard<std::mutex> const lock(ResonantTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        int blocklengths[] = {1, 1, 3, 1, 1};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Aint displacements[] = {offsetof(ResonantTerm, ResCoeff),
                                    offsetof(ResonantTerm, NonResCoeff),
                                    offsetof(ResonantTerm, Poles),
                                    offsetof(ResonantTerm, isz1z2),
                                    offsetof(ResonantTerm, Weight)};
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        MPI_Datatype types[] = {
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType ResCoeff
            POMEROL_MPI_DOUBLE_COMPLEX, // ComplexType NonResCoeff
            MPI_DOUBLE,                 // RealType Poles[3]
            POMEROL_MPI_BOOL,           // bool isz1z2
            MPI_LONG                    // long Weight
        };
        MPI_Type_create_struct(5, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

//
// TwoParticleGFPart
//
TwoParticleGFPart::TwoParticleGFPart(MonomialOperatorPart const& O1,
                                     MonomialOperatorPart const& O2,
                                     MonomialOperatorPart const& O3,
                                     MonomialOperatorPart const& CX4,
                                     HamiltonianPart const& Hpart1,
                                     HamiltonianPart const& Hpart2,
                                     HamiltonianPart const& Hpart3,
                                     HamiltonianPart const& Hpart4,
                                     DensityMatrixPart const& DMpart1,
                                     DensityMatrixPart const& DMpart2,
                                     DensityMatrixPart const& DMpart3,
                                     DensityMatrixPart const& DMpart4,
                                     Permutation3 Permutation)
    : Thermal(DMpart1.beta),
      ComputableObject(),
      O1(O1),
      O2(O2),
      O3(O3),
      CX4(CX4),
      Hpart1(Hpart1),
      Hpart2(Hpart2),
      Hpart3(Hpart3),
      Hpart4(Hpart4),
      DMpart1(DMpart1),
      DMpart2(DMpart2),
      DMpart3(DMpart3),
      DMpart4(DMpart4),
      Permutation(std::move(Permutation)),
      NonResonantTerms(NonResonantTerm::Hash(), NonResonantTerm::KeyEqual(), NonResonantTerm::IsNegligible()),
      ResonantTerms(ResonantTerm::Hash(), ResonantTerm::KeyEqual(), ResonantTerm::IsNegligible()) {}

void TwoParticleGFPart::compute() {
    if(getStatus() >= Computed)
        return;

    if(O1.isComplex() || O2.isComplex() || O3.isComplex() || CX4.isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();
}

template <bool Complex> void TwoParticleGFPart::computeImpl() {
    NonResonantTerms.clear();
    ResonantTerms.clear();

    RealType beta = DMpart1.beta;
    // I don't have any pen now, so I'm writing here:
    // <1 | O1 | 2> <2 | O2 | 3> <3 | O3 |4> <4| CX4 |1>
    // Iterate over all values of |1><1| and |3><3|
    // Chase indices |2> and <2|, |4> and <4|.
    RowMajorMatrixType<Complex> const& O1matrix = O1.getRowMajorValue<Complex>();
    ColMajorMatrixType<Complex> const& O2matrix = O2.getColMajorValue<Complex>();
    RowMajorMatrixType<Complex> const& O3matrix = O3.getRowMajorValue<Complex>();
    ColMajorMatrixType<Complex> const& CX4matrix = CX4.getColMajorValue<Complex>();

    InnerQuantumState index1Max =
        CX4matrix.outerSize(); // One can not make a cutoff in external index for evaluating 2PGF
    InnerQuantumState index3Max = O2matrix.outerSize();

    std::vector<InnerQuantumState> Index4List;
    Index4List.reserve(index1Max * index3Max);

    for(InnerQuantumState index1 = 0; index1 < index1Max; ++index1)
        for(InnerQuantumState index3 = 0; index3 < index3Max; ++index3) {
            typename ColMajorMatrixType<Complex>::InnerIterator index4bra_iter(CX4matrix, index1);
            typename RowMajorMatrixType<Complex>::InnerIterator index4ket_iter(O3matrix, index3);
            Index4List.clear();
            while(index4bra_iter && index4ket_iter) {
                if(chaseIndices<Complex>(index4ket_iter, index4bra_iter)) {
                    Index4List.push_back(index4bra_iter.index());
                    ++index4bra_iter;
                    ++index4ket_iter;
                }
            };

            if(!Index4List.empty()) {
                RealType E1 = Hpart1.getEigenValue(index1);
                RealType E3 = Hpart3.getEigenValue(index3);
                RealType weight1 = DMpart1.getWeight(index1);
                RealType weight3 = DMpart3.getWeight(index3);

                typename ColMajorMatrixType<Complex>::InnerIterator index2bra_iter(O2matrix, index3);
                typename RowMajorMatrixType<Complex>::InnerIterator index2ket_iter(O1matrix, index1);
                while(index2bra_iter && index2ket_iter) {
                    if(chaseIndices<Complex>(index2ket_iter, index2bra_iter)) {

                        InnerQuantumState index2 = index2ket_iter.index();
                        RealType E2 = Hpart2.getEigenValue(index2);
                        RealType weight2 = DMpart2.getWeight(index2);

                        for(unsigned long index4 : Index4List) {
                            RealType E4 = Hpart4.getEigenValue(index4);
                            RealType weight4 = DMpart4.getWeight(index4);
                            if(weight1 + weight2 + weight3 + weight4 >= CoefficientTolerance) {
                                ComplexType MatrixElement = index2ket_iter.value() * index2bra_iter.value() *
                                                            O3matrix.coeff(index3, index4) *
                                                            CX4matrix.coeff(index4, index1);

                                MatrixElement *= Permutation.sign;

                                addMultiterm(MatrixElement, beta, E1, E2, E3, E4, weight1, weight2, weight3, weight4);
                            }
                        }
                        ++index2bra_iter;
                        ++index2ket_iter;
                    }
                }
            }
        }

    INFO("Total " << NonResonantTerms.size() << "+" << ResonantTerms.size() << "="
                  << NonResonantTerms.size() + ResonantTerms.size() << " terms");

    assert(NonResonantTerms.check_terms());
    assert(ResonantTerms.check_terms());

    setStatus(Computed);
}

inline void TwoParticleGFPart::addMultiterm(ComplexType Coeff,
                                            RealType beta,
                                            RealType Ei,
                                            RealType Ej,
                                            RealType Ek,
                                            RealType El,
                                            RealType Wi,
                                            RealType Wj,
                                            RealType Wk,
                                            RealType Wl) {
    RealType P1 = Ej - Ei;
    RealType P2 = Ek - Ej;
    RealType P3 = El - Ek;

    // Non-resonant part of the multiterm
    ComplexType CoeffZ2 = -Coeff * (Wj + Wk);
    if(std::abs(CoeffZ2) > CoefficientTolerance)
        NonResonantTerms.add_term(NonResonantTerm(CoeffZ2, P1, P2, P3, false));
    ComplexType CoeffZ4 = Coeff * (Wi + Wl);
    if(std::abs(CoeffZ4) > CoefficientTolerance)
        NonResonantTerms.add_term(NonResonantTerm(CoeffZ4, P1, P2, P3, true));

    // Resonant part of the multiterm
    ComplexType CoeffZ1Z2Res = Coeff * beta * Wi;
    ComplexType CoeffZ1Z2NonRes = Coeff * (Wk - Wi);
    if(std::abs(CoeffZ1Z2Res) > CoefficientTolerance || abs(CoeffZ1Z2NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(ResonantTerm(CoeffZ1Z2Res, CoeffZ1Z2NonRes, P1, P2, P3, true));
    ComplexType CoeffZ2Z3Res = -Coeff * beta * Wj;
    ComplexType CoeffZ2Z3NonRes = Coeff * (Wj - Wl);
    if(std::abs(CoeffZ2Z3Res) > CoefficientTolerance || abs(CoeffZ2Z3NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(ResonantTerm(CoeffZ2Z3Res, CoeffZ2Z3NonRes, P1, P2, P3, false));
}

ComplexType TwoParticleGFPart::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const {
    long MatsubaraNumberOdd1 = 2 * MatsubaraNumber1 + 1;
    long MatsubaraNumberOdd2 = 2 * MatsubaraNumber2 + 1;
    long MatsubaraNumberOdd3 = 2 * MatsubaraNumber3 + 1;
    return (*this)(MatsubaraSpacing * RealType(MatsubaraNumberOdd1),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd2),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd3));
}

ComplexType TwoParticleGFPart::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const {
    std::array<ComplexType, 3> Frequencies = {z1, z2, -z3};

    z1 = Frequencies[Permutation.perm[0]];
    z2 = Frequencies[Permutation.perm[1]];
    z3 = Frequencies[Permutation.perm[2]];

    if(getStatus() != Computed) {
        throw StatusMismatch(
            "2PGFPart: Calling operator() on uncomputed container. Did you purge all the terms when called compute()?");
    }

    return NonResonantTerms(z1, z2, z3) + ResonantTerms(z1, z2, z3, ReduceResonanceTolerance);
}

void TwoParticleGFPart::clear() {
    NonResonantTerms.clear();
    ResonantTerms.clear();
    setStatus(Constructed);
}

} // namespace Pomerol
