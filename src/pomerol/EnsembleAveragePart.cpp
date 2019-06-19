#include "pomerol/EnsembleAveragePart.h"

namespace Pomerol{

EnsembleAveragePart::EnsembleAveragePart( const QuadraticOperatorPart& A,
                                        const HamiltonianPart& Hpart, const DensityMatrixPart& DMpart) :
                                        Thermal(DMpart),
                                        Hpart(Hpart), DMpart(DMpart),
                                        A(A),
                                        result(0),
                                        MatrixElementTolerance(1e-8)
{}

void EnsembleAveragePart::compute(void)
{
    result = 0;

//    A.print_to_screen();
    // Blocks (submatrices) of A
    const RowMajorMatrixType& Amatrix = A.getRowMajorValue();
    QuantumState outerSize = Amatrix.outerSize();

    // Iterate over all values of the outer index.
//    for(QuantumState index1=0; index1<outerSize; ++index1){
//        // <index1|A|Ainner>
//        RowMajorMatrixType::InnerIterator Ainner(Amatrix,index1);
//
//        // While we are not at the last column of Amatrix or at the last row of CXmatrix.
//        while(Ainner){
//            QuantumState A_index2 = Ainner.index();
//
//            // diagonal element
//            if(A_index2 == index1)
//                result += Ainner.value() * DMpart.getWeight(index1);
//            ++Ainner;   // The next non-zero element
//        }
//    }

    // Sum up <index1|A|index1> * weight(index1)
    for(QuantumState index1=0; index1<outerSize; ++index1)
        result += Amatrix.coeff(index1, index1) * DMpart.getWeight(index1);
}

} // end of namespace Pomerol
