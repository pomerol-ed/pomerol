#include "pomerol/EnsembleAverage.hpp"

namespace Pomerol {

EnsembleAverage::EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                                 const MonomialOperator& A, const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), DM(DM)
{}

EnsembleAverage::EnsembleAverage(const EnsembleAverage& EA) :
    Thermal(EA.beta), ComputableObject(EA), S(EA.S), H(EA.H), A(EA.A), DM(EA.DM), Result(EA.Result)
{}

void EnsembleAverage::compute()
{
    if(getStatus() >= Prepared) return;

    // Find out non-trivial blocks of A.
    MonomialOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(auto Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++) {
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright){
            // check if retained blocks are included. If not, do not push.
            if(DM.isRetained(Aleft)) {
                if(A.isComplex())
                    Result += computeImpl<true>((MonomialOperatorPart&)A.getPartFromLeftIndex(Aleft), DM.getPart(Aleft));
                else
                    Result += computeImpl<false>((MonomialOperatorPart&)A.getPartFromLeftIndex(Aleft), DM.getPart(Aleft));
            }
        }
    }

    setStatus(Prepared);
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
    for(Eigen::Index Index = 0; Index < Amatrix.outerSize(); ++Index)
        result_part += Amatrix.coeff(Index, Index) * DMpart.getWeight(Index);
    return result_part;
}

} // namespace Pomerol
