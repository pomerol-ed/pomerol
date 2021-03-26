#include "pomerol/FieldOperatorPart.h"

using std::stringstream;

namespace Pomerol{

void prune_if_real(RowMajorMatrixType<true> &, RealType) {}
void prune_if_real(RowMajorMatrixType<false> & a, RealType tol) {
  a.prune(tol);
}

template<bool Complex>
FieldOperatorPart<Complex>::FieldOperatorPart(
        const IndexClassification<Complex> &IndexInfo,
        const StatesClassification<Complex> &S,
        const HamiltonianPart<Complex> &HFrom,
        const HamiltonianPart<Complex> &HTo,
        ParticleIndex PIndex) :
        ComputableObject(), IndexInfo(IndexInfo), S(S), HFrom(HFrom), HTo(HTo), PIndex(PIndex),
        MatrixElementTolerance(1e-8)
{}

template<bool Complex>
void FieldOperatorPart<Complex>::compute()
{
    if ( Status >= Computed ) return;
    BlockNumber to = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();

    const std::vector<FockState>& toStates = S.getFockStates(to);
    const std::vector<FockState>& fromStates = S.getFockStates(from);

    MatrixType<Complex> RightMat(fromStates.size(), fromStates.size());
    MatrixType<Complex> LeftMat(toStates.size(), fromStates.size());
    RightMat.setZero();
    LeftMat.setZero();

    /* Rotation is done in the following way:
     * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
     * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
     * We use the fact each column of O_{lk} has only one nonzero elements.
     * */
    for (std::vector<FockState>::const_iterator CurrentState = fromStates.begin();
                                                CurrentState < fromStates.end(); CurrentState++) {
	    FockState K=*CurrentState;
        std::map<FockState, MelemType<Complex>> result1 = O->actRight(K);
        if (result1.size()) {
            FockState L=result1.begin()->first;
            int sign = int(real(result1.begin()->second));
	        if ( L!=ERROR_FOCK_STATE && std::abs(sign)>std::numeric_limits<RealType>::epsilon() ) {
		        InnerQuantumState l=S.getInnerState(L), k=S.getInnerState(K);

                for (InnerQuantumState n=0; n<toStates.size(); n++) {
                    LeftMat(n,k) = conj(HTo.getMatrixElement(l,n));
                }

                for (InnerQuantumState m=0; m<fromStates.size(); m++) {
                    RightMat(k,m) = RealType(sign) * HFrom.getMatrixElement(k,m);
                }
	        }
        }
    }

// Workaround for Eigen issue 1224
// https://gitlab.com/libeigen/eigen/-/issues/1224
//
// Affected versions are some betas of 3.3 but not the 3.3 release
#if EIGEN_VERSION_AT_LEAST(3,2,90) && EIGEN_MAJOR_VERSION<3
    elementsRowMajor = MatrixType(LeftMat * RightMat).sparseView(MatrixElementTolerance);
#else
    elementsRowMajor = (LeftMat * RightMat).sparseView(MatrixElementTolerance);
#endif

    prune_if_real(elementsRowMajor, MatrixElementTolerance);
    elementsColMajor = elementsRowMajor;
    Status = Computed;
}

template<bool Complex>
const ColMajorMatrixType<Complex>& FieldOperatorPart<Complex>::getColMajorValue(void) const
{
    return elementsColMajor;
}

template<bool Complex>
const RowMajorMatrixType<Complex>& FieldOperatorPart<Complex>::getRowMajorValue(void) const
{
    return elementsRowMajor;
}

template<bool Complex>
void FieldOperatorPart<Complex>::print_to_screen() const  //print to screen C and CX
{
    BlockNumber to   = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();
    INFO(S.getQuantumNumbers(from) << "->" << S.getQuantumNumbers(to));
    for (size_t P=0; P<elementsColMajor.outerSize(); ++P)
    for(typename ColMajorMatrixType<Complex>::InnerIterator it(elementsColMajor,P); it; ++it) {
	    FockState N = S.getFockState(to, it.row());
	    FockState M = S.getFockState(from, it.col());
	    INFO(N <<" " << M << " : " << it.value());
    }
}

template<bool Complex>
BlockNumber FieldOperatorPart<Complex>::getLeftIndex(void) const
{
    return HTo.getBlockNumber();
}

template<bool Complex>
BlockNumber FieldOperatorPart<Complex>::getRightIndex(void) const
{
    return HFrom.getBlockNumber();
}

// Specialized methods

template<bool Complex>
AnnihilationOperatorPart<Complex>::AnnihilationOperatorPart(const IndexClassification<Complex> &IndexInfo,
                                                            const StatesClassification<Complex> &S,
                                                            const HamiltonianPart<Complex> &HFrom,
                                                            const HamiltonianPart<Complex> &HTo,
                                                            ParticleIndex PIndex) :
    FieldOperatorPart<Complex>(IndexInfo,S,HFrom,HTo,PIndex)
{
    this->O = new Pomerol::OperatorPresets::C<Complex>(PIndex);
}

template<bool Complex>
CreationOperatorPart<Complex>::CreationOperatorPart(const IndexClassification<Complex> &IndexInfo,
                                                    const StatesClassification<Complex> &S,
                                                    const HamiltonianPart<Complex> &HFrom,
                                                    const HamiltonianPart<Complex> &HTo,
                                                    ParticleIndex PIndex) :
    FieldOperatorPart<Complex>(IndexInfo,S,HFrom,HTo,PIndex)
{
    this->O = new Pomerol::OperatorPresets::Cdag<Complex>(PIndex);
}

template<bool Complex>
const CreationOperatorPart<Complex>& AnnihilationOperatorPart<Complex>::transpose() const
{
    CreationOperatorPart<Complex> *CX = new CreationOperatorPart<Complex>(this->IndexInfo,
                                                                          this->S,
                                                                          this->HTo,
                                                                          this->HFrom,
                                                                          this->PIndex); // swapped h_to and h_from
    CX->elementsRowMajor = this->elementsRowMajor.transpose();
    CX->elementsColMajor = this->elementsColMajor.transpose();
    return *CX;
}

template<bool Complex>
const AnnihilationOperatorPart<Complex>& CreationOperatorPart<Complex>::transpose() const
{
    AnnihilationOperatorPart<Complex> *C = new AnnihilationOperatorPart<Complex>(this->IndexInfo,
                                                                                 this->S,
                                                                                 this->HTo,
                                                                                 this->HFrom,
                                                                                 this->PIndex); // swapped h_to and h_from
    C->elementsRowMajor = this->elementsRowMajor.transpose();
    C->elementsColMajor = this->elementsColMajor.transpose();
    return *C;
}

template<bool Complex>
QuadraticOperatorPart<Complex>::QuadraticOperatorPart(const IndexClassification<Complex> &IndexInfo,
                                                      const StatesClassification<Complex> &S,
                                                      const HamiltonianPart<Complex> &HFrom,
                                                      const HamiltonianPart<Complex> &HTo,
                                                      ParticleIndex PIndex1, ParticleIndex PIndex2) :
        FieldOperatorPart<Complex>(IndexInfo,S,HFrom,HTo,9999), Index1(PIndex1), Index2(PIndex2)
        // Index=9999 dummy
{
    this->O = new Pomerol::OperatorPresets::N_offdiag<Complex>(PIndex1, PIndex2);
}

// Explicit instantiations: Real case

template class FieldOperatorPart<false>;
template class AnnihilationOperatorPart<false>;
template class CreationOperatorPart<false>;
template class QuadraticOperatorPart<false>;

// Explicit instantiations: Complex case

template class FieldOperatorPart<true>;
template class AnnihilationOperatorPart<true>;
template class CreationOperatorPart<true>;
template class QuadraticOperatorPart<true>;

} // end of namespace Pomerol
