#include "pomerol/LibcommuteEigen.h"
#include "pomerol/MonomialOperatorPart.h"

#include <libcommute/loperator/mapped_basis_view.hpp>

namespace Pomerol {

void MonomialOperatorPart::compute()
{
    if(Status >= Computed) return;
    if(isComplex() && HFrom.isComplex())
        computeImpl<true, true>();
    else if(isComplex() && !HFrom.isComplex())
        computeImpl<true, false>();
    else
        computeImpl<false, false>();
    Status = Computed;
}

template<bool Complex, bool HComplex>
void MonomialOperatorPart::computeImpl() {
    BlockNumber to = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();

    const std::vector<QuantumState>& toStates = S.getFockStates(to);
    const std::vector<QuantumState>& fromStates = S.getFockStates(from);

    /* Rotation is done in the following way:
    * O_{nm} = \sum_{lk} U^{+}_{nl} O_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
    * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
    * We use the fact each column of O_{lk} has only one nonzero elements.
    * */
    MatrixType<Complex> OURight(toStates.size(), fromStates.size());

    auto fromMapper = libcommute::basis_mapper(fromStates);
    auto toMapper = libcommute::basis_mapper(toStates);

    auto const& U = HFrom.getMatrix<HComplex>();

    auto const& MOp_ = *static_cast<const LOperatorType<Complex>*>(MOp);

    for(InnerQuantumState st = 0; st < fromStates.size(); ++st) {
        auto fromView = fromMapper.make_const_view_no_ref(U.col(st));
        auto toView = toMapper.make_view_no_ref(OURight.col(st));
        MOp_(fromView, toView);
    }

    auto const& ULeft = HTo.getMatrix<HComplex>().adjoint();

// Workaround for Eigen issue 1224
// https://gitlab.com/libeigen/eigen/-/issues/1224
//
// Affected versions are some betas of 3.3 but not the 3.3 release
#if EIGEN_VERSION_AT_LEAST(3,2,90) && EIGEN_MAJOR_VERSION<3
    elementsRowMajor = std::make_shared<RowMajorMatrixType<Complex>>(
        MatrixType<Complex>(ULeft * OURight).sparseView(MatrixElementTolerance)
    );
#else
    elementsRowMajor = std::make_shared<RowMajorMatrixType<Complex>>(
        (ULeft * OURight).sparseView(MatrixElementTolerance)
    );
#endif

    elementsColMajor = std::make_shared<ColMajorMatrixType<Complex>>(
        *std::static_pointer_cast<const RowMajorMatrixType<Complex>>(elementsRowMajor)
    );
}

void MonomialOperatorPart::setFromAdjoint(const MonomialOperatorPart &part) {
    assert(isComplex() == part.isComplex());
    if(Status >= Computed) return;

    if(isComplex()) {
        elementsRowMajor = std::make_shared<RowMajorMatrixType<true>>(part.getColMajorValue<true>().adjoint());
        elementsColMajor = std::make_shared<RowMajorMatrixType<true>>(part.getRowMajorValue<true>().adjoint());
    } else {
        elementsRowMajor = std::make_shared<RowMajorMatrixType<false>>(part.getColMajorValue<false>().adjoint());
        elementsColMajor = std::make_shared<RowMajorMatrixType<false>>(part.getRowMajorValue<false>().adjoint());
    }

    Status = Computed;
}

template<bool Complex>
ColMajorMatrixType<Complex>& MonomialOperatorPart::getColMajorValue()
{
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<ColMajorMatrixType<Complex>>(elementsColMajor);
}
template ColMajorMatrixType<true>& MonomialOperatorPart::getColMajorValue<true>();
template ColMajorMatrixType<false>& MonomialOperatorPart::getColMajorValue<false>();

template<bool Complex>
const ColMajorMatrixType<Complex>& MonomialOperatorPart::getColMajorValue() const
{
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const ColMajorMatrixType<Complex>>(elementsColMajor);
}
template const ColMajorMatrixType<true>& MonomialOperatorPart::getColMajorValue<true>() const;
template const ColMajorMatrixType<false>& MonomialOperatorPart::getColMajorValue<false>() const;

template<bool Complex>
RowMajorMatrixType<Complex>& MonomialOperatorPart::getRowMajorValue()
{
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<RowMajorMatrixType<Complex>>(elementsRowMajor);
}
template RowMajorMatrixType<true>& MonomialOperatorPart::getRowMajorValue<true>();
template RowMajorMatrixType<false>& MonomialOperatorPart::getRowMajorValue<false>();

template<bool Complex>
const RowMajorMatrixType<Complex>& MonomialOperatorPart::getRowMajorValue() const
{
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const RowMajorMatrixType<Complex>>(elementsRowMajor);
}
template const RowMajorMatrixType<true>& MonomialOperatorPart::getRowMajorValue<true>() const;
template const RowMajorMatrixType<false>& MonomialOperatorPart::getRowMajorValue<false>() const;

void MonomialOperatorPart::print_to_screen() const
{
    if(Complex)
        print_to_screenImpl<true>();
    else
        print_to_screenImpl<false>();
}

template<bool Complex>
void MonomialOperatorPart::print_to_screenImpl() const {
    BlockNumber to   = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();
    auto const& mat = getColMajorValue<Complex>();

    for(size_t P = 0; P < mat.outerSize(); ++P)
    for(typename ColMajorMatrixType<Complex>::InnerIterator it(mat, P); it; ++it) {
        QuantumState N = S.getFockState(to, it.row());
        QuantumState M = S.getFockState(from, it.col());
        INFO(N <<" " << M << " : " << it.value());
    }
}

BlockNumber MonomialOperatorPart::getLeftIndex() const
{
    return HTo.getBlockNumber();
}

BlockNumber MonomialOperatorPart::getRightIndex() const
{
    return HFrom.getBlockNumber();
}

} // end of namespace Pomerol
