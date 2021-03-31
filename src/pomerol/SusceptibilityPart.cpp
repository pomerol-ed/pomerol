#include "pomerol/SusceptibilityPart.h"

namespace Pomerol{

template<bool Complex>
SusceptibilityPart<Complex>::Term::Term(ComplexType Residue, RealType Pole) :
    Residue(Residue), Pole(Pole) {};

// BOSON: minus sign before Residue according to the definition \chi(tau) = < A(tau) B >
template<bool Complex>
ComplexType SusceptibilityPart<Complex>::Term::operator()(ComplexType Frequency) const { return -Residue/(Frequency - Pole); }

template<bool Complex>
ComplexType SusceptibilityPart<Complex>::Term::operator()(RealType tau, RealType beta) const {
    return Pole > 0 ? Residue*exp(-tau*Pole)/(1 - exp(-beta*Pole)) :
                      Residue*exp((beta-tau)*Pole)/(exp(beta*Pole) - 1);
}

template<bool Complex>
inline auto SusceptibilityPart<Complex>::Term::operator+=(const Term& AnotherTerm) -> Term&
{
    Residue += AnotherTerm.Residue;
    return *this;
}

template<bool Complex>
std::ostream& operator<<(std::ostream& out, typename SusceptibilityPart<Complex>::Term const& T)
{
    out << T.Residue << "/(z - " << T.Pole << ")";
    return out;
}

template<bool Complex>
SusceptibilityPart<Complex>::SusceptibilityPart( const QuadraticOperatorPart<Complex>& A, const QuadraticOperatorPart<Complex>& B,
                                        const HamiltonianPart<Complex>& HpartInner, const HamiltonianPart<Complex>& HpartOuter,
                                        const DensityMatrixPart<Complex>& DMpartInner, const DensityMatrixPart<Complex>& DMpartOuter) :
                                        Thermal(DMpartInner),
                                        Terms(typename Term::Compare(1e-8), typename Term::IsNegligible(1e-8)),
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        A(A), B(B),
                                        MatrixElementTolerance(1e-8),
                                        ReduceResonanceTolerance(1e-8),
                                        ReduceTolerance(1e-8),
                                        ZeroPoleWeight(0)
{}

template<bool Complex>
void SusceptibilityPart<Complex>::compute(void)
{
    Terms.clear();

    // Blocks (submatrices) of A and B
    RowMajorMatrixType<Complex> const& Amatrix = A.getRowMajorValue();
    ColMajorMatrixType<Complex> const& Bmatrix = B.getColMajorValue();
    QuantumState outerSize = Amatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Amatrix and empty columns of Bmatrix.
    for(QuantumState index1=0; index1<outerSize; ++index1){
        // <index1|A|Ainner><Binner|B|index1>
        typename RowMajorMatrixType<Complex>::InnerIterator Ainner(Amatrix,index1);
        typename ColMajorMatrixType<Complex>::InnerIterator Binner(Bmatrix,index1);

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

template class SusceptibilityPart<false>;
template class SusceptibilityPart<true>;

} // end of namespace Pomerol
