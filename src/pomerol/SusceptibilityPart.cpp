#include "pomerol/SusceptibilityPart.hpp"

#include <cassert>

namespace Pomerol {

SusceptibilityPart::Term::Term(ComplexType Residue, RealType Pole) :
    Residue(Residue), Pole(Pole) {}

// BOSON: minus sign before Residue according to the definition \chi(tau) = < A(tau) B >
ComplexType SusceptibilityPart::Term::operator()(ComplexType Frequency) const
{
    return -Residue/(Frequency - Pole);
}

ComplexType SusceptibilityPart::Term::operator()(RealType tau, RealType beta) const {
    using std::exp;
    return Pole > 0 ? Residue*exp(-tau*Pole)/(1 - exp(-beta*Pole)) :
                      Residue*exp((beta-tau)*Pole)/(exp(beta*Pole) - 1);
}

inline
SusceptibilityPart::Term& SusceptibilityPart::Term::operator+=(const Term& AnotherTerm)
{
    Residue += AnotherTerm.Residue;
    return *this;
}

SusceptibilityPart::SusceptibilityPart( const MonomialOperatorPart& A, const MonomialOperatorPart& B,
                                        const HamiltonianPart& HpartInner, const HamiltonianPart& HpartOuter,
                                        const DensityMatrixPart& DMpartInner, const DensityMatrixPart& DMpartOuter) :
                                        Thermal(DMpartInner),
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        A(A), B(B),
                                        Terms(Term::Compare(), Term::IsNegligible())
{}

void SusceptibilityPart::compute()
{
    if(A.isComplex()) {
        if(B.isComplex())
            computeImpl<true, true>();
        else
            computeImpl<true, false>();
    } else {
        if(B.isComplex())
            computeImpl<false, true>();
        else
            computeImpl<false, false>();
    }
}

template<bool AComplex, bool BComplex>
void SusceptibilityPart::computeImpl()
{
    Terms.clear();

    // Blocks (submatrices) of A and B
    const RowMajorMatrixType<AComplex>& Amatrix = A.getRowMajorValue<AComplex>();
    const ColMajorMatrixType<BComplex>& Bmatrix = B.getColMajorValue<BComplex>();
    QuantumState outerSize = Amatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Amatrix and empty columns of Bmatrix.
    for(QuantumState index1 = 0; index1 < outerSize; ++index1) {
        // <index1|A|Ainner><Binner|B|index1>
        typename RowMajorMatrixType<AComplex>::InnerIterator Ainner(Amatrix,index1);
        typename ColMajorMatrixType<BComplex>::InnerIterator Binner(Bmatrix,index1);

        // While we are not at the last column of Amatrix or at the last row of Bmatrix.
        while(Ainner && Binner){
            QuantumState A_index2 = Ainner.index();
            QuantumState B_index2 = Binner.index();

            // A meaningful matrix element
            if(A_index2 == B_index2){
                RealType Pole = HpartInner.getEigenValue(A_index2) - HpartOuter.getEigenValue(index1);
                if(std::abs(Pole) < ReduceResonanceTolerance){
                    // BOSON: pole at zero energy
                    ZeroPoleWeight += Ainner.value() * Binner.value() * DMpartOuter.getWeight(index1);
                }
                else{
                    // BOSON: minus sign before the second term
                    ComplexType Residue = Ainner.value() * Binner.value() *
                                          (DMpartOuter.getWeight(index1) - DMpartInner.getWeight(A_index2));
                    if(std::abs(Residue) > MatrixElementTolerance) // Is the residue relevant?
                    {
                        // Create a new term and append it to the list.
                        Terms.add_term(Term(Residue, Pole));
                    }
                }
                ++Ainner; // The next non-zero element
                ++Binner; // The next non-zero element
            } else {
                // Chasing: one index runs down the other index
                if(B_index2 < A_index2) for(;QuantumState(Binner.index())<A_index2; ++Binner);
                else for(;QuantumState(Ainner.index())<B_index2; ++Ainner);
            }
        }
    }

    assert(Terms.check_terms());
}

} // namespace Pomerol
