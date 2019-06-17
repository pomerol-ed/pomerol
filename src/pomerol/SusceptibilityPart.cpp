#include "pomerol/SusceptibilityPart.h"

namespace Pomerol{

SusceptibilityPart::Term::Term(ComplexType Residue, RealType Pole) :
    Residue(Residue), Pole(Pole) {};
ComplexType SusceptibilityPart::Term::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

ComplexType SusceptibilityPart::Term::operator()(RealType tau, RealType beta) const {
    return Pole > 0 ? -Residue*exp(-tau*Pole)/(1 + exp(-beta*Pole)) :
                      -Residue*exp((beta-tau)*Pole)/(exp(beta*Pole) + 1);
}

inline
SusceptibilityPart::Term& SusceptibilityPart::Term::operator+=(const Term& AnotherTerm)
{
    Residue += AnotherTerm.Residue;
    return *this;
}

std::ostream& operator<<(std::ostream& out, const SusceptibilityPart::Term& T)
{
    out << T.Residue << "/(z - " << T.Pole << ")";
    return out;
}

SusceptibilityPart::SusceptibilityPart( const QuadraticOperatorPart& A, const QuadraticOperatorPart& B,
                                        const HamiltonianPart& HpartInner, const HamiltonianPart& HpartOuter,
                                        const DensityMatrixPart& DMpartInner, const DensityMatrixPart& DMpartOuter) :
                                        Thermal(DMpartInner),
                                        Terms(Term::Compare(1e-8), Term::IsNegligible(1e-8)),
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        A(A), B(B),
                                        MatrixElementTolerance(1e-8),
                                        ReduceResonanceTolerance(1e-8),
                                        ReduceTolerance(1e-8)
{}

void SusceptibilityPart::compute(void)
{
    Terms.clear();

    // Blocks (submatrices) of A and B
    const RowMajorMatrixType& Amatrix = A.getRowMajorValue();
    const ColMajorMatrixType& Bmatrix = B.getColMajorValue();
    QuantumState outerSize = Amatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Amatrix and empty columns of Bmatrix.
    for(QuantumState index1=0; index1<outerSize; ++index1){
        // <index1|A|Ainner><Binner|B|index1>
        RowMajorMatrixType::InnerIterator Ainner(Amatrix,index1);
        ColMajorMatrixType::InnerIterator Binner(Bmatrix,index1);

        // While we are not at the last column of Amatrix or at the last row of Bmatrix.
        while(Ainner && Binner){
            QuantumState A_index2 = Ainner.index();
            QuantumState B_index2 = Binner.index();

            // A meaningful matrix element
            if(A_index2 == B_index2){
                ComplexType Residue = Ainner.value() * Binner.value() *
                                      (DMpartOuter.getWeight(index1) + DMpartInner.getWeight(A_index2));
                if(abs(Residue) > MatrixElementTolerance) // Is the residue relevant?
                {
                    // Create a new term and append it to the list.
                    RealType Pole = HpartInner.getEigenValue(A_index2) - HpartOuter.getEigenValue(index1);
                    Terms.add_term(Term(Residue, Pole));
                    //DEBUG("<" << C.S.getFockState(HpartInner.getBlockNumber(), A_index2) << "|" << Residue << "|" <<  C.S.getFockState(HpartInner.getBlockNumber(),B_index2) << ">" );
                };
                ++Ainner;   // The next non-zero element
                ++Binner;  // The next non-zero element
            }else{
                // Chasing: one index runs down the other index
                if(B_index2 < A_index2) for(;QuantumState(Binner.index())<A_index2; ++Binner);
                else for(;QuantumState(Ainner.index())<B_index2; ++Ainner);
            }
        }
    }

    assert(Terms.check_terms());
}

} // end of namespace Pomerol
