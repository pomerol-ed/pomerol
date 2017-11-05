#include "pomerol/FieldOperatorPart.h"

using std::stringstream;

namespace Pomerol{

FieldOperatorPart::FieldOperatorPart(
        const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom,  const HamiltonianPart &HTo, ParticleIndex PIndex) :
        ComputableObject(), IndexInfo(IndexInfo), S(S), HFrom(HFrom), HTo(HTo), PIndex(PIndex),
        MatrixElementTolerance(1e-8)
{}

void FieldOperatorPart::compute()
{
    if ( Status >= Computed ) return;
    BlockNumber to = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();

    const std::vector<FockState>& toStates = S.getFockStates(to);
    const std::vector<FockState>& fromStates = S.getFockStates(from);
    MatrixType tempElements(toStates.size(), fromStates.size());
    tempElements.setZero();

    /* Rotation is done in the following way:
     * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
     * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates. */
    for (std::vector<FockState>::const_iterator CurrentState = fromStates.begin();
                                                CurrentState < fromStates.end(); CurrentState++) {
	    FockState K=*CurrentState;
        std::map<FockState, MelemType> result1 = O->actRight(K);
        if (result1.size()) {
            FockState L=result1.begin()->first;
            #ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
            int sign = int(std::real(result1.begin()->second));
            #else
            int sign = result1.begin()->second;
            #endif
	        if ( L!=ERROR_FOCK_STATE && std::abs(sign)>std::numeric_limits<RealType>::epsilon() ) {
		        InnerQuantumState l=S.getInnerState(L), k=S.getInnerState(K);

                MatrixType Un(toStates.size(),1);
                for (InnerQuantumState n=0; n<toStates.size(); n++) {
#ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
                    Un(n,0) = std::conj(HTo.getMatrixElement(l,n));
#else
                    Un(n,0) = HTo.getMatrixElement(l,n);
#endif
                }

                MatrixType Um(1,fromStates.size());
                for (InnerQuantumState m=0; m<fromStates.size(); m++) {
                    Um(0,m) = HFrom.getMatrixElement(k,m);
                }
                tempElements += RealType(sign) * Un * Um;
	        }
        }
    }

//    for (std::vector<QuantumState>::const_iterator CurrentState = toStates.begin();
//                                                   CurrentState < toStates.end(); CurrentState++)
    elementsRowMajor = RowMajorMatrixType(toStates.size(),fromStates.size());
    elementsRowMajor = tempElements.sparseView(MatrixElementTolerance);
    #ifndef POMEROL_COMPLEX_MATRIX_ELEMENTS
    elementsRowMajor.prune(MatrixElementTolerance);
    #endif
    elementsColMajor = elementsRowMajor;
    Status = Computed;
}

const ColMajorMatrixType& FieldOperatorPart::getColMajorValue(void) const
{
    return elementsColMajor;
}

const RowMajorMatrixType& FieldOperatorPart::getRowMajorValue(void) const
{
    return elementsRowMajor;
}

void FieldOperatorPart::print_to_screen() const  //print to screen C and CX
{
    BlockNumber to   = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();
    INFO(S.getQuantumNumbers(from) << "->" << S.getQuantumNumbers(to));
    for (size_t P=0; P<elementsColMajor.outerSize(); ++P)
	for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it) {
	    FockState N = S.getFockState(to, it.row());
	    FockState M = S.getFockState(from, it.col());
	    INFO(N <<" " << M << " : " << it.value());
        }
}

BlockNumber FieldOperatorPart::getLeftIndex(void) const
{
    return HTo.getBlockNumber();
}

BlockNumber FieldOperatorPart::getRightIndex(void) const
{
    return HFrom.getBlockNumber();
}

// Specialized methods

AnnihilationOperatorPart::AnnihilationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S,
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{
    O = new Pomerol::OperatorPresets::C(PIndex);
}

CreationOperatorPart::CreationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S,
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{
    O = new Pomerol::OperatorPresets::Cdag(PIndex);
}

const CreationOperatorPart& AnnihilationOperatorPart::transpose() const
{
    CreationOperatorPart *CX = new CreationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    CX->elementsRowMajor = elementsRowMajor.transpose();
    CX->elementsColMajor = elementsColMajor.transpose();
    return *CX;
}

const AnnihilationOperatorPart& CreationOperatorPart::transpose() const
{
    AnnihilationOperatorPart *C = new AnnihilationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    C->elementsRowMajor = elementsRowMajor.transpose();
    C->elementsColMajor = elementsColMajor.transpose();
    return *C;
}

} // end of namespace Pomerol
