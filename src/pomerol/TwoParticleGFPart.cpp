//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


#include "pomerol/TwoParticleGFPart.h"

namespace Pomerol{

// Make the lagging index catch up or outrun the leading index.
inline bool chaseIndices(RowMajorMatrixType::InnerIterator& index1_iter,
                         ColMajorMatrixType::InnerIterator& index2_iter)
{
    InnerQuantumState index1 = index1_iter.index();
    InnerQuantumState index2 = index2_iter.index();

    if(index1 == index2) return true;

    if(index1 < index2)
        for(;InnerQuantumState(index1_iter.index())<index2 && index1_iter; ++index1_iter);
    else
        for(;InnerQuantumState(index2_iter.index())<index1 && index2_iter; ++index2_iter);

    return false;
}

//
// TwoParticleGFPart::NonResonantTerm
//
inline
TwoParticleGFPart::NonResonantTerm::NonResonantTerm(ComplexType Coeff, RealType P1, RealType P2, RealType P3, bool isz4) :
Coeff(Coeff), isz4(isz4)
{
    Poles[0] = P1; Poles[1] = P2; Poles[2] = P3; Weight=1;
}

inline
TwoParticleGFPart::NonResonantTerm& TwoParticleGFPart::NonResonantTerm::operator+=(
                    const NonResonantTerm& AnotherTerm)
{
    long combinedWeight=Weight + AnotherTerm.Weight;
    for (unsigned short p=0; p<3; ++p) Poles[p]= (Weight*Poles[p] + AnotherTerm.Weight*AnotherTerm.Poles[p])/combinedWeight;
    Weight=combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

//
// TwoParticleGFPart::ResonantTerm
//
inline
TwoParticleGFPart::ResonantTerm::ResonantTerm(ComplexType ResCoeff, ComplexType NonResCoeff,
                                              RealType P1, RealType P2, RealType P3, bool isz1z2):
ResCoeff(ResCoeff), NonResCoeff(NonResCoeff), isz1z2(isz1z2)
{
    Poles[0] = P1; Poles[1] = P2; Poles[2] = P3; Weight=1;
}

inline
TwoParticleGFPart::ResonantTerm& TwoParticleGFPart::ResonantTerm::operator+=(
                const ResonantTerm& AnotherTerm)
{
    long combinedWeight=Weight + AnotherTerm.Weight;
    for (unsigned short p=0; p<3; ++p) Poles[p]= (Weight*Poles[p] + AnotherTerm.Weight*AnotherTerm.Poles[p])/combinedWeight;
    Weight=combinedWeight;
    ResCoeff += AnotherTerm.ResCoeff;
    NonResCoeff += AnotherTerm.NonResCoeff;
    return *this;
}

//
// TwoParticleGFPart
//
TwoParticleGFPart::TwoParticleGFPart(
                const FieldOperatorPart& O1, const FieldOperatorPart& O2,
                const FieldOperatorPart& O3, const CreationOperatorPart& CX4,
                const HamiltonianPart& Hpart1, const HamiltonianPart& Hpart2,
                const HamiltonianPart& Hpart3, const HamiltonianPart& Hpart4,
                const DensityMatrixPart& DMpart1, const DensityMatrixPart& DMpart2,
                const DensityMatrixPart& DMpart3, const DensityMatrixPart& DMpart4,
                Permutation3 Permutation) :
    Thermal(DMpart1),
    ComputableObject(),
    NonResonantTerms(NonResonantTerm::Compare(1e-8), NonResonantTerm::IsNegligible(1e-16)),
    ResonantTerms(ResonantTerm::Compare(1e-8), ResonantTerm::IsNegligible(1e-16)),
    O1(O1), O2(O2), O3(O3), CX4(CX4),
    Hpart1(Hpart1), Hpart2(Hpart2), Hpart3(Hpart3), Hpart4(Hpart4),
    DMpart1(DMpart1), DMpart2(DMpart2), DMpart3(DMpart3), DMpart4(DMpart4),
    Permutation(Permutation),
    KroneckerSymbolTolerance(1e-16),
    ReduceResonanceTolerance(1e-8),
    CoefficientTolerance (1e-16),
    MultiTermCoefficientTolerance (1e-5)
{}

void TwoParticleGFPart::compute()
{
    NonResonantTerms.clear();
    ResonantTerms.clear();

    RealType beta = DMpart1.beta;
    // I don't have any pen now, so I'm writing here:
    // <1 | O1 | 2> <2 | O2 | 3> <3 | O3 |4> <4| CX4 |1>
    // Iterate over all values of |1><1| and |3><3|
    // Chase indices |2> and <2|, |4> and <4|.
    const RowMajorMatrixType& O1matrix = O1.getRowMajorValue();
    const ColMajorMatrixType& O2matrix = O2.getColMajorValue();
    const RowMajorMatrixType& O3matrix = O3.getRowMajorValue();
    const ColMajorMatrixType& CX4matrix = CX4.getColMajorValue();

    InnerQuantumState index1;
    InnerQuantumState index1Max = CX4matrix.outerSize(); // One can not make a cutoff in external index for evaluating 2PGF

    InnerQuantumState index3;
    InnerQuantumState index3Max = O2matrix.outerSize();

    std::vector<InnerQuantumState> Index4List;
    Index4List.reserve(index1Max*index3Max);

    for(index1=0; index1<index1Max; ++index1)
    for(index3=0; index3<index3Max; ++index3){
        ColMajorMatrixType::InnerIterator index4bra_iter(CX4matrix,index1);
        RowMajorMatrixType::InnerIterator index4ket_iter(O3matrix,index3);
        Index4List.clear();
        while (index4bra_iter && index4ket_iter){
            if(chaseIndices(index4ket_iter,index4bra_iter)){
                Index4List.push_back(index4bra_iter.index());
                ++index4bra_iter;
                ++index4ket_iter;
            }
        };

        if (!Index4List.empty())
        {
            RealType E1 = Hpart1.getEigenValue(index1);
            RealType E3 = Hpart3.getEigenValue(index3);
            RealType weight1 = DMpart1.getWeight(index1);
            RealType weight3 = DMpart3.getWeight(index3);

            ColMajorMatrixType::InnerIterator index2bra_iter(O2matrix,index3);
            RowMajorMatrixType::InnerIterator index2ket_iter(O1matrix,index1);
            while (index2bra_iter && index2ket_iter){
                if (chaseIndices(index2ket_iter,index2bra_iter)){

                    InnerQuantumState index2 = index2ket_iter.index();
                    RealType E2 = Hpart2.getEigenValue(index2);
                    RealType weight2 = DMpart2.getWeight(index2);

                    for (unsigned long p4 = 0; p4 < Index4List.size(); ++p4)
                    {
                        InnerQuantumState index4 = Index4List[p4];//*pIndex4;
                        RealType E4 = Hpart4.getEigenValue(index4);
                        RealType weight4 = DMpart4.getWeight(index4);
                        if (weight1 + weight2 + weight3 + weight4 >= CoefficientTolerance) {
                            ComplexType MatrixElement = index2ket_iter.value()*
                                                        index2bra_iter.value()*
                                                        O3matrix.coeff(index3,index4)*
                                                        CX4matrix.coeff(index4,index1);

                            MatrixElement *= Permutation.sign;

                            addMultiterm(MatrixElement,beta,E1,E2,E3,E4,weight1,weight2,weight3,weight4);
                        }
                    }
                    ++index2bra_iter;
                    ++index2ket_iter;
                };
            }
        };
    }

    std::cout << "Total " << NonResonantTerms.size() << "+" << ResonantTerms.size() << "="
              << NonResonantTerms.size() + ResonantTerms.size() << " terms" << std::endl << std::flush;

    assert(NonResonantTerms.check_terms());
    assert(ResonantTerms.check_terms());

    Status = Computed;
}

inline
void TwoParticleGFPart::addMultiterm(ComplexType Coeff, RealType beta,
                      RealType Ei, RealType Ej, RealType Ek, RealType El,
                      RealType Wi, RealType Wj, RealType Wk, RealType Wl)
{
    RealType P1 = Ej - Ei;
    RealType P2 = Ek - Ej;
    RealType P3 = El - Ek;

    // Non-resonant part of the multiterm
    ComplexType CoeffZ2 = -Coeff*(Wj + Wk);
    if(abs(CoeffZ2) > CoefficientTolerance)
        NonResonantTerms.add_term(
            NonResonantTerm(CoeffZ2,P1,P2,P3,false));
    ComplexType CoeffZ4 = Coeff*(Wi + Wl);
    if(abs(CoeffZ4) > CoefficientTolerance)
        NonResonantTerms.add_term(
            NonResonantTerm(CoeffZ4,P1,P2,P3,true));

    // Resonant part of the multiterm
    ComplexType CoeffZ1Z2Res = Coeff*beta*Wi;
    ComplexType CoeffZ1Z2NonRes = Coeff*(Wk - Wi);
    if(abs(CoeffZ1Z2Res) > CoefficientTolerance || abs(CoeffZ1Z2NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(
            ResonantTerm(CoeffZ1Z2Res,CoeffZ1Z2NonRes,P1,P2,P3,true));
    ComplexType CoeffZ2Z3Res = -Coeff*beta*Wj;
    ComplexType CoeffZ2Z3NonRes = Coeff*(Wj - Wl);
    if(abs(CoeffZ2Z3Res) > CoefficientTolerance || abs(CoeffZ2Z3NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(
            ResonantTerm(CoeffZ2Z3Res,CoeffZ2Z3NonRes,P1,P2,P3,false));
}

size_t TwoParticleGFPart::getNumNonResonantTerms() const
{
    return NonResonantTerms.size();
}

size_t TwoParticleGFPart::getNumResonantTerms() const
{
    return ResonantTerms.size();
}

const Permutation3& TwoParticleGFPart::getPermutation() const
{
    return Permutation;
}

ComplexType TwoParticleGFPart::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long MatsubaraNumberOdd1 = 2*MatsubaraNumber1 + 1;
    long MatsubaraNumberOdd2 = 2*MatsubaraNumber2 + 1;
    long MatsubaraNumberOdd3 = 2*MatsubaraNumber3 + 1;
    return (*this)(MatsubaraSpacing * RealType(MatsubaraNumberOdd1),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd2),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd3));
}

ComplexType TwoParticleGFPart::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const
{
    ComplexType Frequencies[3] = {  z1, z2, -z3 };

    z1 = Frequencies[Permutation.perm[0]];
    z2 = Frequencies[Permutation.perm[1]];
    z3 = Frequencies[Permutation.perm[2]];

    if (Status != Computed) {
        throw std::logic_error("2PGFPart : Calling operator() on empty container, did you purge all the terms when called compute()");
    }

    return NonResonantTerms(z1, z2, z3) + ResonantTerms(z1, z2, z3, KroneckerSymbolTolerance);
}

const TermList<TwoParticleGFPart::ResonantTerm>& TwoParticleGFPart::getResonantTerms() const
{
    return ResonantTerms;
}

const TermList<TwoParticleGFPart::NonResonantTerm>& TwoParticleGFPart::getNonResonantTerms() const
{
    return NonResonantTerms;
}

void TwoParticleGFPart::clear()
{
    NonResonantTerms.clear();
    ResonantTerms.clear();
    Status = Constructed;
}

} // end of namespace Pomerol

