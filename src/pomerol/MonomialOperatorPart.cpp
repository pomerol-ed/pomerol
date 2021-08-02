#include "pomerol/MonomialOperatorPart.hpp"

#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/loperator/mapped_basis_view.hpp>

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Pomerol {

void MonomialOperatorPart::compute()
{
    if(getStatus() >= Computed) return;

    if(MOpComplex && HFrom.isComplex())
        computeImpl<true, true>();
    else if(MOpComplex && !HFrom.isComplex())
        computeImpl<true, false>();
    else if(!MOpComplex && HFrom.isComplex())
        computeImpl<false, true>();
    else
        computeImpl<false, false>();

    setStatus(Computed);
}

template<bool MOpC, bool HC>
void MonomialOperatorPart::computeImpl() {
    constexpr bool C = MOpC || HC;

    BlockNumber to = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();

    const std::vector<QuantumState>& toStates = S.getFockStates(to);
    const std::vector<QuantumState>& fromStates = S.getFockStates(from);

    /* Rotation is done in the following way:
    * O_{nm} = \sum_{lk} U^{+}_{nl} O_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
    * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
    * We use the fact each column of O_{lk} has only one nonzero elements.
    * */
    MatrixType<C> OURight(toStates.size(), fromStates.size());

    auto fromMapper = libcommute::basis_mapper(fromStates);
    auto toMapper = libcommute::basis_mapper(toStates);

    auto const& U = HFrom.getMatrix<HC>();

    auto const& MOp_ = *static_cast<const LOperatorTypeRC<MOpC>*>(MOp);

    for(InnerQuantumState st = 0; st < fromStates.size(); ++st) {
        auto fromView = fromMapper.make_const_view_no_ref(U.col(st));
        auto toView = toMapper.make_view_no_ref(OURight.col(st));
        MOp_(fromView, toView);
    }

    auto const& ULeft = HTo.getMatrix<HC>().adjoint();

// Workaround for Eigen issue 1224
// https://gitlab.com/libeigen/eigen/-/issues/1224
//
// Affected versions are some betas of 3.3 but not the 3.3 release
#if EIGEN_VERSION_AT_LEAST(3,2,90) && EIGEN_MAJOR_VERSION<3
    elementsRowMajor = std::make_shared<RowMajorMatrixType<C>>(
        MatrixType<C>(ULeft * OURight).sparseView(MatrixElementTolerance)
    );
#else
    elementsRowMajor = std::make_shared<RowMajorMatrixType<C>>(
        (ULeft * OURight).sparseView(MatrixElementTolerance)
    );
#endif

    elementsColMajor = std::make_shared<ColMajorMatrixType<C>>(
        *std::static_pointer_cast<const RowMajorMatrixType<C>>(elementsRowMajor)
    );
}

void MonomialOperatorPart::setFromAdjoint(const MonomialOperatorPart &part) {
    assert(isComplex() == part.isComplex());
    assert(getLeftIndex() == part.getRightIndex());
    assert(getRightIndex() == part.getLeftIndex());

    if(getStatus() >= Computed) return;

    if(isComplex()) {
        elementsRowMajor = std::make_shared<RowMajorMatrixType<true>>(part.getColMajorValue<true>().adjoint());
        elementsColMajor = std::make_shared<ColMajorMatrixType<true>>(part.getRowMajorValue<true>().adjoint());
    } else {
        elementsRowMajor = std::make_shared<RowMajorMatrixType<false>>(part.getColMajorValue<false>().adjoint());
        elementsColMajor = std::make_shared<ColMajorMatrixType<false>>(part.getRowMajorValue<false>().adjoint());
    }

    setStatus(Computed);
}

template<bool C>
ColMajorMatrixType<C>& MonomialOperatorPart::getColMajorValue()
{
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<ColMajorMatrixType<C>>(elementsColMajor);
}
template ColMajorMatrixType<true>& MonomialOperatorPart::getColMajorValue<true>();
template ColMajorMatrixType<false>& MonomialOperatorPart::getColMajorValue<false>();

template<bool C>
const ColMajorMatrixType<C>& MonomialOperatorPart::getColMajorValue() const
{
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const ColMajorMatrixType<C>>(elementsColMajor);
}
template const ColMajorMatrixType<true>& MonomialOperatorPart::getColMajorValue<true>() const;
template const ColMajorMatrixType<false>& MonomialOperatorPart::getColMajorValue<false>() const;

template<bool C>
RowMajorMatrixType<C>& MonomialOperatorPart::getRowMajorValue()
{
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<RowMajorMatrixType<C>>(elementsRowMajor);
}
template RowMajorMatrixType<true>& MonomialOperatorPart::getRowMajorValue<true>();
template RowMajorMatrixType<false>& MonomialOperatorPart::getRowMajorValue<false>();

template<bool C>
const RowMajorMatrixType<C>& MonomialOperatorPart::getRowMajorValue() const
{
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const RowMajorMatrixType<C>>(elementsRowMajor);
}
template const RowMajorMatrixType<true>& MonomialOperatorPart::getRowMajorValue<true>() const;
template const RowMajorMatrixType<false>& MonomialOperatorPart::getRowMajorValue<false>() const;

template<bool C>
void MonomialOperatorPart::streamOutputImpl(std::ostream & os) const {
    BlockNumber to   = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();
    auto const& mat = getColMajorValue<C>();

    for(size_t P = 0; P < mat.outerSize(); ++P) {
        for(typename ColMajorMatrixType<C>::InnerIterator it(mat, P); it; ++it) {
            QuantumState N = S.getFockState(to, it.row());
            QuantumState M = S.getFockState(from, it.col());
            os << N << " " << M << " : " << it.value() << std::endl;
        }
    }
}
template void MonomialOperatorPart::streamOutputImpl<true>(std::ostream & os) const;
template void MonomialOperatorPart::streamOutputImpl<false>(std::ostream & os) const;

} // namespace Pomerol
