#include "pomerol/GreensFunctionPart.h"

namespace Pomerol{

template<bool Complex>
GreensFunctionPart<Complex>::Term::Term(ComplexType Residue, RealType Pole) :
    Residue(Residue), Pole(Pole) {};

template<bool Complex>
ComplexType GreensFunctionPart<Complex>::Term::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

template<bool Complex>
ComplexType GreensFunctionPart<Complex>::Term::operator()(RealType tau, RealType beta) const {
    return Pole > 0 ? -Residue*exp(-tau*Pole)/(1 + exp(-beta*Pole)) :
                      -Residue*exp((beta-tau)*Pole)/(exp(beta*Pole) + 1);
}

template<bool Complex>
inline auto GreensFunctionPart<Complex>::Term::operator+=(const Term& AnotherTerm) -> Term&
{
    Residue += AnotherTerm.Residue;
    return *this;
}

template<bool Complex>
GreensFunctionPart<Complex>::GreensFunctionPart(const AnnihilationOperatorPart<Complex>& C,
                                                const CreationOperatorPart<Complex>& CX,
                                                const HamiltonianPart<Complex>& HpartInner,
                                                const HamiltonianPart<Complex>& HpartOuter,
                                                const DensityMatrixPart<Complex>& DMpartInner,
                                                const DensityMatrixPart<Complex>& DMpartOuter) :
                                        Thermal(DMpartInner),
                                        Terms(typename Term::Compare(1e-8), typename Term::IsNegligible(1e-8)),
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        C(C), CX(CX),
                                        MatrixElementTolerance(1e-8),
                                        ReduceResonanceTolerance(1e-8),
                                        ReduceTolerance(1e-8)
{}

template<bool Complex>
void GreensFunctionPart<Complex>::compute(void)
{
    Terms.clear();

    // Blocks (submatrices) of C and CX
    const RowMajorMatrixType<Complex>& Cmatrix = C.getRowMajorValue();
    const ColMajorMatrixType<Complex>& CXmatrix = CX.getColMajorValue();
    QuantumState outerSize = Cmatrix.outerSize();

    // Iterate over all values of the outer index.
    // TODO: should be optimized - skip empty rows of Cmatrix and empty columns of CXmatrix.
    for(QuantumState index1=0; index1<outerSize; ++index1){
        // <index1|C|Cinner><CXinner|CX|index1>
        typename RowMajorMatrixType<Complex>::InnerIterator Cinner(Cmatrix,index1);
        typename ColMajorMatrixType<Complex>::InnerIterator CXinner(CXmatrix,index1);

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
                    Terms.add_term(Term(Residue, Pole));
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

    assert(Terms.check_terms());
}

template class GreensFunctionPart<false>;
template class GreensFunctionPart<true>;

} // end of namespace Pomerol
