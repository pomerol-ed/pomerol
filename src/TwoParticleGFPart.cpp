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


#include "TwoParticleGFPart.h"

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

inline
bool TwoParticleGFPart::NonResonantTerm::isSimilarTo(const NonResonantTerm& AnotherTerm) const
{   
    return isz4 == AnotherTerm.isz4 && 
           (fabs(Poles[0] - AnotherTerm.Poles[0]) < ReduceResonanceTolerance) &&
           (fabs(Poles[1] - AnotherTerm.Poles[1]) < ReduceResonanceTolerance) &&
           (fabs(Poles[2] - AnotherTerm.Poles[2]) < ReduceResonanceTolerance);
}

//
// TwoParticleGFPart::ResonantTerm
//
inline
TwoParticleGFPart::ResonantTerm::ResonantTerm(ComplexType ResCoeff, ComplexType NonResCoeff,
                                              RealType P1, RealType P2, RealType P3, bool isz1z2) : 
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

inline
bool TwoParticleGFPart::ResonantTerm::isSimilarTo(const ResonantTerm& AnotherTerm) const
{   
    return isz1z2 == AnotherTerm.isz1z2 &&
           (fabs(Poles[0] - AnotherTerm.Poles[0]) < ReduceResonanceTolerance) &&
           (fabs(Poles[1] - AnotherTerm.Poles[1]) < ReduceResonanceTolerance) &&
           (fabs(Poles[2] - AnotherTerm.Poles[2]) < ReduceResonanceTolerance);
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
    O1(O1), O2(O2), O3(O3), CX4(CX4), 
    Hpart1(Hpart1), Hpart2(Hpart2), Hpart3(Hpart3), Hpart4(Hpart4),
    DMpart1(DMpart1), DMpart2(DMpart2), DMpart3(DMpart3), DMpart4(DMpart4),
    Permutation(Permutation)
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
    InnerQuantumState index1Max = CX4matrix.outerSize();

    InnerQuantumState index3;
    InnerQuantumState index3Max = O2matrix.outerSize();

    std::list<InnerQuantumState> Index4List;

    unsigned long ResonantTermsUnreducedSize=0;
    unsigned long NonResonantTermsUnreducedSize=0;
    unsigned long ResonantTermsPreviousSize=0;
    unsigned long NonResonantTermsPreviousSize=0;

    for(index1=0; index1<index1Max; ++index1){
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

                    for (std::list<InnerQuantumState>::iterator pIndex4 = Index4List.begin(); pIndex4!=Index4List.end(); ++pIndex4) 
                    {
                        InnerQuantumState index4 = *pIndex4;
                        RealType E4 = Hpart4.getEigenValue(index4);                       
                        RealType weight4 = DMpart4.getWeight(index4);

                        ComplexType MatrixElement = index2ket_iter.value()*
                                                    index2bra_iter.value()*
                                                    O3matrix.coeff(index3,index4)*
                                                    CX4matrix.coeff(index4,index1);

                        MatrixElement *= Permutation.sign;

                        addMultiterm(MatrixElement,beta,E1,E2,E3,E4,weight1,weight2,weight3,weight4);
                    }
                    ++index2bra_iter;
                    ++index2ket_iter;
                };
            }
        };
    }
    if (ResonantTerms.size()-ResonantTermsPreviousSize + NonResonantTerms.size() - NonResonantTermsPreviousSize > ReduceInvocationThreshold ){ 
        INFO_NONEWLINE(NonResonantTerms.size()-NonResonantTermsPreviousSize << " nonresonant + " << ResonantTerms.size() - ResonantTermsPreviousSize << " resonant = ");
        INFO_NONEWLINE(ResonantTerms.size()-ResonantTermsPreviousSize + NonResonantTerms.size() - NonResonantTermsPreviousSize);
        INFO_NONEWLINE(" terms reduced to ");

        NonResonantTermsUnreducedSize+=NonResonantTerms.size();
        ResonantTermsUnreducedSize+= ResonantTerms.size();

        RealType NonResonantTolerance = MultiTermCoefficientTolerance*(index1+1)/NonResonantTermsUnreducedSize/(index1Max+1);
        RealType ResonantTolerance = MultiTermCoefficientTolerance*(index1+1)/ResonantTermsUnreducedSize/(index1Max+1);

        this->reduceTerms(NonResonantTolerance, ResonantTolerance, NonResonantTerms, ResonantTerms); 
        NonResonantTermsPreviousSize = NonResonantTerms.size();
        ResonantTermsPreviousSize = ResonantTerms.size(); 

        INFO_NONEWLINE(NonResonantTermsPreviousSize << "+" << ResonantTermsPreviousSize << " = ");
        INFO(NonResonantTermsPreviousSize + ResonantTermsPreviousSize << " with tolerances: " << NonResonantTolerance << ", " << ResonantTolerance);
        };
    };

    NonResonantTermsUnreducedSize=(NonResonantTermsUnreducedSize>0)?NonResonantTermsUnreducedSize:NonResonantTerms.size();
    ResonantTermsUnreducedSize=(ResonantTermsUnreducedSize>0)?ResonantTermsUnreducedSize:ResonantTerms.size();
    if (ResonantTermsUnreducedSize + NonResonantTermsUnreducedSize > 0){
        INFO_NONEWLINE("Total " << NonResonantTermsUnreducedSize << " nonresonant + " << ResonantTermsUnreducedSize << " resonant = ");
        INFO_NONEWLINE(NonResonantTermsUnreducedSize+ResonantTermsUnreducedSize << " terms reduced to ");
        this->reduceTerms(MultiTermCoefficientTolerance/(NonResonantTermsUnreducedSize+1), MultiTermCoefficientTolerance/(ResonantTermsUnreducedSize+1), NonResonantTerms, ResonantTerms);
        INFO_NONEWLINE(NonResonantTerms.size() << "+" << ResonantTerms.size() << " = ");
        INFO(NonResonantTerms.size() + ResonantTerms.size()  << " with tolerances: " << MultiTermCoefficientTolerance/(NonResonantTermsUnreducedSize+1) << ", " << MultiTermCoefficientTolerance/(ResonantTermsUnreducedSize+1));

    };
}


void TwoParticleGFPart::reduceTerms(const RealType NonResonantTolerance, const RealType ResonantTolerance, std::list<NonResonantTerm> &NonResonantTerms, std::list<ResonantTerm>& ResonantTerms)
{
    #ifndef noReduction
    // Sieve reduction of the non-resonant terms
    //int i=0
    for(std::list<NonResonantTerm>::iterator it1 = NonResonantTerms.begin(); it1 != NonResonantTerms.end();){
        std::list<NonResonantTerm>::iterator it2 = it1;
    //    DEBUG(++i << "/" << NonResonantTerms.size());
        for(it2++; it2 != NonResonantTerms.end();){
            if(it1->isSimilarTo(*it2)){
                *it1 += *it2;
                it2 = NonResonantTerms.erase(it2);
            }else
                it2++;
        }
        
        if(abs(it1->Coeff) < NonResonantTolerance)
            it1 = NonResonantTerms.erase(it1);
        else
            it1++;
    }
    
    // Sieve reduction of the resonant terms
    for(std::list<ResonantTerm>::iterator it1 = ResonantTerms.begin(); it1 != ResonantTerms.end();){
        std::list<ResonantTerm>::iterator it2 = it1;
        for(it2++; it2 != ResonantTerms.end();){
            if(it1->isSimilarTo(*it2)){
                *it1 += *it2;
                it2 = ResonantTerms.erase(it2);
            }else
                it2++;
        }
        
        if(abs(it1->ResCoeff) + abs(it1->NonResCoeff) < ResonantTolerance)
            it1 = ResonantTerms.erase(it1);
        else
            it1++;
    }
    #else 
        INFO("No reduction is done due to -DnoReduction flag turned on");
    #endif // endif :: #ifndef noReduction
    //DEBUG("After: " << NonResonantTerms.size() << " non resonant terms + " << ResonantTerms.size() << " resonant terms = " << ResonantTerms.size()+NonResonantTerms.size());
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
        NonResonantTerms.push_back(
            NonResonantTerm(CoeffZ2,P1,P2,P3,false));
    ComplexType CoeffZ4 = Coeff*(Wi + Wl);
    if(abs(CoeffZ4) > CoefficientTolerance)
        NonResonantTerms.push_back(
            NonResonantTerm(CoeffZ4,P1,P2,P3,true));

    // Resonant part of the multiterm
    ComplexType CoeffZ1Z2Res = Coeff*beta*Wi;
    ComplexType CoeffZ1Z2NonRes = Coeff*(Wk - Wi);
    if(abs(CoeffZ1Z2Res) > CoefficientTolerance || abs(CoeffZ1Z2NonRes) > CoefficientTolerance)
        ResonantTerms.push_back(
            ResonantTerm(CoeffZ1Z2Res,CoeffZ1Z2NonRes,P1,P2,P3,true));   
    ComplexType CoeffZ2Z3Res = -Coeff*beta*Wj;
    ComplexType CoeffZ2Z3NonRes = Coeff*(Wj - Wl);
    if(abs(CoeffZ2Z3Res) > CoefficientTolerance || abs(CoeffZ2Z3NonRes) > CoefficientTolerance)
        ResonantTerms.push_back(
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
    ComplexType Frequencies[3] = {  MatsubaraSpacing * RealType(MatsubaraNumberOdd1),
                                    MatsubaraSpacing * RealType(MatsubaraNumberOdd2),
                                   -MatsubaraSpacing * RealType(MatsubaraNumberOdd3)};

    ComplexType z1 = Frequencies[Permutation.perm[0]];
    ComplexType z2 = Frequencies[Permutation.perm[1]];
    ComplexType z3 = Frequencies[Permutation.perm[2]];

    ComplexType Value = 0;
    for(std::list<NonResonantTerm>::const_iterator pTerm = NonResonantTerms.begin(); pTerm != NonResonantTerms.end(); ++pTerm)
        Value += (*pTerm)(z1,z2,z3);
    for(std::list<ResonantTerm>::const_iterator pTerm = ResonantTerms.begin(); pTerm != ResonantTerms.end(); ++pTerm)
        Value += (*pTerm)(z1,z2,z3);

    return Value;
}

const std::list<TwoParticleGFPart::NonResonantTerm>& TwoParticleGFPart::getNonResonantTerms() const
{
    return NonResonantTerms;
}

const std::list<TwoParticleGFPart::ResonantTerm>& TwoParticleGFPart::getResonantTerms() const
{
    return ResonantTerms;
}

void TwoParticleGFPart::clear()
{
    NonResonantTerms.clear();
    ResonantTerms.clear();
}

} // end of namespace Pomerol

