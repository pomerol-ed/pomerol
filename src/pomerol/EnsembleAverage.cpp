#include "pomerol/EnsembleAverage.h"

namespace Pomerol{

EnsembleAverage::EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                                 const MonomialOperator& A, const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), DM(DM), result(0)
{}

EnsembleAverage::EnsembleAverage(const EnsembleAverage& EA) :
    Thermal(EA.beta), ComputableObject(EA), S(EA.S), H(EA.H), A(EA.A), DM(EA.DM), result(EA.result)
{}

void EnsembleAverage::compute()
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A.
    MonomialOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(MonomialOperator::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++){
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright){
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) ){
                if(A.isComplex())
                    result += computeImpl<true>((MonomialOperatorPart&)A.getPartFromLeftIndex(Aleft), DM.getPart(Aleft));
                else
                    result += computeImpl<false>((MonomialOperatorPart&)A.getPartFromLeftIndex(Aleft), DM.getPart(Aleft));
            }
        }
    }

    Status = Prepared;
}

// This function is called directly in prepare()
template<bool Complex>
ComplexType EnsembleAverage::computeImpl(const MonomialOperatorPart& Apart,
                                         const DensityMatrixPart& DMpart)
{
    // Blocks (submatrices) of A
    const RowMajorMatrixType<Complex>& Amatrix = Apart.getRowMajorValue<Complex>();

    // Sum up <index1|A|index1> * weight(index1)
    ComplexType result_part = 0;
    for(QuantumState index1=0; index1<Amatrix.outerSize(); ++index1)
        result_part += Amatrix.coeff(index1, index1) * DMpart.getWeight(index1);
    return result_part;
}

} // end of namespace Pomerol
