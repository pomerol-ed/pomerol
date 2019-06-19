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

    // Blocks (submatrices) of A
    const RowMajorMatrixType& Amatrix = A.getRowMajorValue();

    // Sum up <index1|A|index1> * weight(index1)
    for(QuantumState index1=0; index1<Amatrix.outerSize(); ++index1)
        result += Amatrix.coeff(index1, index1) * DMpart.getWeight(index1);
}

} // end of namespace Pomerol
