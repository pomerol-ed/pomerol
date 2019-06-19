#include "pomerol/EnsembleAverage.h"

namespace Pomerol{

EnsembleAverage::EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                                 const QuadraticOperator& A, const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), DM(DM), result(0)
{}

EnsembleAverage::EnsembleAverage(const EnsembleAverage& EA) :
    Thermal(EA.beta), ComputableObject(EA), S(EA.S), H(EA.H), A(EA.A), DM(EA.DM), result(EA.result)
{}

void EnsembleAverage::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A.
    FieldOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(FieldOperator::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++){
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright){
            DEBUG(S.getQuantumNumbers(Aleft) << "|" << S.getQuantumNumbers(Aright) );
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) ){
                result += compute((QuadraticOperatorPart&)A.getPartFromLeftIndex(Aleft),
                                  H.getPart(Aleft), DM.getPart(Aleft));
//                EnsembleAveragePart part((QuadraticOperatorPart&)A.getPartFromLeftIndex(Aleft),
//                                         H.getPart(Aleft), DM.getPart(Aleft));
//                part.compute();
//                result += part.getResult();
            }
        }
    }

    Status = Prepared;
}

// This function is called directly in prepare()
ComplexType EnsembleAverage::compute(const QuadraticOperatorPart& Apart,
                                     const HamiltonianPart& Hpart,
                                     const DensityMatrixPart& DMpart)
{
    // Blocks (submatrices) of A
    const RowMajorMatrixType& Amatrix = Apart.getRowMajorValue();

    // Sum up <index1|A|index1> * weight(index1)
    ComplexType result_part = 0;
    for(QuantumState index1=0; index1<Amatrix.outerSize(); ++index1)
        result_part += Amatrix.coeff(index1, index1) * DMpart.getWeight(index1);
    return result_part;
}

} // end of namespace Pomerol
