/** \file src/GreensFunctionPart.cpp
** \brief Part of a Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "GreensFunctionPart.h"

GreensFunctionPart::GreensTerm::GreensTerm(ComplexType Residue, ComplexType Pole) : Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::GreensTerm::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

std::ostream& operator<<(std::ostream& out, const GreensFunctionPart::GreensTerm& Term)
{
    out << Term.Residue << "/(z - " << Term.Pole << ")";
    return out;
}

GreensFunctionPart::GreensFunctionPart( AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                                        HamiltonianPart& HpartInner, HamiltonianPart& HpartOuter,
                                        DensityMatrixPart& DMpartInner, DensityMatrixPart& DMpartOuter) :
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        C(C), CX(CX)
{}

void GreensFunctionPart::compute(void)
{
    Terms.clear();

    // Blocks (submatrices) of C and CX
    RowMajorMatrixType& Cmatrix = C.getRowMajorValue();
    ColMajorMatrixType& CXmatrix = CX.getColMajorValue();
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
                                      (DMpartOuter.weight(index1) + DMpartInner.weight(C_index2));
                if(abs(Residue) > MatrixElementTolerance) // Is the residue relevant?
                {
                    // Create a new term and append it to the list.
                    ComplexType Pole = HpartInner.reV(C_index2) - HpartOuter.reV(index1);
                    Terms.push_back(GreensTerm(Residue,Pole));
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
}

ComplexType GreensFunctionPart::operator()(long MatsubaraNumber) const
{
    // TODO: Place this variable to a wider scope?
    ComplexType MatsubaraSpacing = I*M_PI/DMpartInner.getBeta();

    ComplexType G = 0;
    for(std::list<GreensTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
        G += (*pTerm)(MatsubaraSpacing*RealType(2*MatsubaraNumber+1));

    return G;
}
