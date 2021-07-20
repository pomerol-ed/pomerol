#include "pomerol/HamiltonianPart.hpp"

#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/loperator/mapped_basis_view.hpp>

#include <Eigen/Eigenvalues>

#include <limits>
#include <sstream>
#include <stdexcept>

//
// class HamiltonianPart
//

namespace Pomerol {

template<bool C>
void HamiltonianPart::initHMatrix() {
    InnerQuantumState BlockSize = S.getBlockSize(Block);
    HMatrix = std::make_shared<MatrixType<C>>(BlockSize, BlockSize);
}

void HamiltonianPart::prepare()
{
    if(getStatus() >= Prepared) return;

    if(isComplex())
        prepareImpl<true>();
    else
        prepareImpl<false>();

    setStatus(Prepared);
}

template<bool C> void HamiltonianPart::prepareImpl()
{
    initHMatrix<C>();

    auto const& HOp_ = *static_cast<const LOperatorTypeRC<C>*>(HOp);
    auto & HMatrix_ = getMatrix<C>();

    auto mapper = libcommute::basis_mapper(S.getFockStates(Block));

    auto BlockSize = S.getBlockSize(Block);
    VectorType<C> ket = VectorType<C>::Zero(BlockSize);
    auto ket_view = mapper.make_const_view(ket);

    for(InnerQuantumState st = 0; st < BlockSize; ++st) {
        auto bra_view = mapper.make_view_no_ref(HMatrix_.col(st));
        ket(st) = 1.0;
        HOp_(ket_view, bra_view);
        ket(st) = .0;
    }

    assert((HMatrix_.adjoint() - HMatrix_).array().abs().maxCoeff()
           < 100*std::numeric_limits<RealType>::epsilon());
}

void HamiltonianPart::compute()
{
    if(getStatus() >= Computed) return;

    if(isComplex())
        computeImpl<true>();
    else
        computeImpl<false>();

    setStatus(Computed);
}

template<bool C> void HamiltonianPart::computeImpl()
{
    auto & HMatrix_ = getMatrix<C>();
    if (HMatrix_.rows() == 1) {
        assert (std::abs(HMatrix_(0,0) - std::real(HMatrix_(0,0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << std::real(HMatrix_(0,0));
        HMatrix_(0,0) = 1;
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixType<C>> Solver(HMatrix_, Eigen::ComputeEigenvectors);
        HMatrix_ = Solver.eigenvectors();
        Eigenvalues = Solver.eigenvalues(); // eigenvectors are ready
    }
}

template<bool C> const MatrixType<C>& HamiltonianPart::getMatrix() const {
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const MatrixType<C>>(HMatrix);
}
template const MatrixType<true>& HamiltonianPart::getMatrix<true>() const;
template const MatrixType<false>& HamiltonianPart::getMatrix<false>() const;

template<bool C> MatrixType<C>& HamiltonianPart::getMatrix() {
    if(C != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<MatrixType<C>>(HMatrix);
}
template MatrixType<true>& HamiltonianPart::getMatrix<true>();
template MatrixType<false>& HamiltonianPart::getMatrix<false>();

RealType HamiltonianPart::getEigenValue(InnerQuantumState state) const
{
    if(getStatus() < Computed) throw exStatusMismatch();
    return Eigenvalues(state);
}

const RealVectorType& HamiltonianPart::getEigenValues() const
{
    if(getStatus() < Computed) throw exStatusMismatch();
    return Eigenvalues;
}

InnerQuantumState HamiltonianPart::getSize() const
{
    return S.getBlockSize(Block);
}

template<bool C>
VectorType<C> HamiltonianPart::getEigenState(InnerQuantumState state) const
{
    if(getStatus() < Computed) throw exStatusMismatch();
    return getMatrix<C>()->col(state);
}

RealType HamiltonianPart::getMinimumEigenvalue() const
{
    if(getStatus() < Computed) throw exStatusMismatch();
    return Eigenvalues.minCoeff();
}

bool HamiltonianPart::reduce(RealType ActualCutoff)
{
    if(getStatus() < Computed) throw exStatusMismatch();

    InnerQuantumState counter = 0;
    for (counter=0; counter < (unsigned int)Eigenvalues.size() && Eigenvalues[counter]<=ActualCutoff; ++counter){};
    INFO("Left " << counter << " eigenvalues : ");

    if (counter) {
        INFO(Eigenvalues.head(counter) << std::endl << "_________");
        Eigenvalues = Eigenvalues.head(counter);
        if(isComplex()) {
            auto & HMatrix_ = getMatrix<true>();
            HMatrix_ = HMatrix_.topLeftCorner(counter, counter);
        } else {
            auto & HMatrix_ = getMatrix<false>();
            HMatrix_ = HMatrix_.topLeftCorner(counter, counter);
        }
        return true;
    } else
        return false;
}

} // namespace Pomerol

