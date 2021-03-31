#include "pomerol/EnsembleAverage.h"

namespace Pomerol{

template<bool Complex>
EnsembleAverage<Complex>::EnsembleAverage(const StatesClassification<Complex>& S,
                                          const Hamiltonian<Complex>& H,
                                          const QuadraticOperator<Complex>& A,
                                          const DensityMatrix<Complex>& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), DM(DM), result(0)
{}

template<bool Complex>
EnsembleAverage<Complex>::EnsembleAverage(const EnsembleAverage& EA) :
    Thermal(EA.beta), ComputableObject(EA), S(EA.S), H(EA.H), A(EA.A), DM(EA.DM), result(EA.result)
{}

template<bool Complex>
void EnsembleAverage<Complex>::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A.
    typename FieldOperator<Complex>::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(auto Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++){
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright){
            DEBUG(S.getQuantumNumbers(Aleft) << "|" << S.getQuantumNumbers(Aright) );
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) ){
                result += compute((QuadraticOperatorPart<Complex>&)A.getPartFromLeftIndex(Aleft),
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
template<bool Complex>
ComplexType EnsembleAverage<Complex>::compute(const QuadraticOperatorPart<Complex>& Apart,
                                     const HamiltonianPart<Complex>& Hpart,
                                     const DensityMatrixPart<Complex>& DMpart)
{
    // Blocks (submatrices) of A
    const RowMajorMatrixType<Complex>& Amatrix = Apart.getRowMajorValue();

    // Sum up <index1|A|index1> * weight(index1)
    ComplexType result_part = 0;
    for(QuantumState index1=0; index1<Amatrix.outerSize(); ++index1)
        result_part += Amatrix.coeff(index1, index1) * DMpart.getWeight(index1);
    return result_part;
}

template class EnsembleAverage<false>;
template class EnsembleAverage<true>;

} // end of namespace Pomerol
