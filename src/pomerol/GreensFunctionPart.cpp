//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file src/GreensFunctionPart.cpp
** \brief Part of a Green's function.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/
#include "GreensFunctionPart.h"

namespace Pomerol{

GreensFunctionPart::Term::Term(ComplexType Residue, RealType Pole) : 
    Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::Term::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

ComplexType GreensFunctionPart::Term::of_tau(RealType tau, RealType beta) const {
    return -Residue*exp(-tau*Pole)/(1 + exp(-beta*Pole));
}

inline
GreensFunctionPart::Term& GreensFunctionPart::Term::operator+=(const Term& AnotherTerm)
{
    Residue += AnotherTerm.Residue;
    return *this;
}

inline
bool GreensFunctionPart::Term::isSimilarTo(const Term& AnotherTerm, RealType ReduceResonanceTolerance) const
{
    return (fabs(Pole - AnotherTerm.Pole) < ReduceResonanceTolerance);
}

std::ostream& operator<<(std::ostream& out, const GreensFunctionPart::Term& T)
{
    out << T.Residue << "/(z - " << T.Pole << ")";
    return out;
}

GreensFunctionPart::GreensFunctionPart( const AnnihilationOperatorPart& C, const CreationOperatorPart& CX, 
                                        const HamiltonianPart& HpartInner, const HamiltonianPart& HpartOuter,
                                        const DensityMatrixPart& DMpartInner, const DensityMatrixPart& DMpartOuter) :
                                        Thermal(DMpartInner),
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        C(C), CX(CX),
                                        MatrixElementTolerance(1e-8),
                                        ReduceResonanceTolerance(1e-8),
                                        ReduceTolerance(1e-8) 
{}

void GreensFunctionPart::compute(void)
{
    Terms.clear();

    // Blocks (submatrices) of C and CX
    const RowMajorMatrixType& Cmatrix = C.getRowMajorValue();
    const ColMajorMatrixType& CXmatrix = CX.getColMajorValue();
    QuantumState outerSize = Cmatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Cmatrix and empty columns of CXmatrix.
    for(QuantumState index1=0; index1<outerSize; ++index1){
        // <index1|C|Cinner><CXinner|CX|index1>
        RowMajorMatrixType::InnerIterator Cinner(Cmatrix,index1);
        ColMajorMatrixType::InnerIterator CXinner(CXmatrix,index1);

        // While we are not at the last column of Cmatrix or at the last row of CXmatrix.
        while(Cinner && CXinner){
            QuantumState C_index2 = Cinner.index();
            QuantumState CX_index2 = CXinner.index();

            // A meaningful matrix element
            if(C_index2 == CX_index2){
                ComplexType Residue = Cinner.value() * CXinner.value() * 
                                      (DMpartOuter.getWeight(index1) + DMpartInner.getWeight(C_index2));
                if(abs(Residue) > MatrixElementTolerance) // Is the residue relevant?
                {
                    // Create a new term and append it to the list.
                    RealType Pole = HpartInner.getEigenValue(C_index2) - HpartOuter.getEigenValue(index1);
                    Terms.push_back(Term(Residue,Pole));
                    //DEBUG("<" << C.S.getFockState(HpartInner.getBlockNumber(), C_index2) << "|" << Residue << "|" <<  C.S.getFockState(HpartInner.getBlockNumber(),CX_index2) << ">" );
                };
                ++Cinner;   // The next non-zero element
                ++CXinner;  // The next non-zero element
            }else{
                // Chasing: one index runs down the other index
                if(CX_index2 < C_index2) for(;QuantumState(CXinner.index())<C_index2; ++CXinner);
                else for(;QuantumState(Cinner.index())<CX_index2; ++Cinner);
            }
        }
    }

    reduceTerms(ReduceTolerance/Terms.size(),Terms);
}

void GreensFunctionPart::reduceTerms(const RealType Tolerance, std::list<Term> &Terms)
{
    // Sieve reduction of the terms
    for(std::list<Term>::iterator it1 = Terms.begin(); it1 != Terms.end();){
        std::list<Term>::iterator it2 = it1;
        for(it2++; it2 != Terms.end();){
            if(it1->isSimilarTo(*it2, ReduceResonanceTolerance)){
                *it1 += *it2;
                it2 = Terms.erase(it2);
            }else
                it2++;
        }

        if(abs(it1->Residue) < Tolerance)
            it1 = Terms.erase(it1);
        else
            it1++;
    }
}

} // end of namespace Pomerol
