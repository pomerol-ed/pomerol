#include "pomerol/LibcommuteEigen.h"
#include "pomerol/HamiltonianPart.h"

#include <libcommute/loperator/mapped_basis_view.hpp>

#include <Eigen/Eigenvalues>

#include <limits>
#include <sstream>
#include <stdexcept>

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#endif

//
// class HamiltonianPart
//

namespace Pomerol {

template<bool Complex>
void HamiltonianPart::initHMatrix() {
    InnerQuantumState BlockSize = S.getBlockSize(Block);
    HMatrix = std::make_shared<MatrixType<Complex>>(BlockSize, BlockSize);
}

void HamiltonianPart::prepare()
{
    if(Status >= Prepared) return;
    if(Complex)
        prepareImpl<true>();
    else
        prepareImpl<false>();
    Status = Prepared;
}

template<bool Complex> void HamiltonianPart::prepareImpl()
{
    initHMatrix<Complex>();

    auto const& HOp_ = *static_cast<const LOperatorType<Complex>*>(HOp);
    auto & HMatrix_ = *std::static_pointer_cast<MatrixType<Complex>>(HMatrix);

    auto mapper = libcommute::basis_mapper(S.getFockStates(Block));

    auto BlockSize = S.getBlockSize(Block);
    VectorType<Complex> ket = VectorType<Complex>::Zero(BlockSize);
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
    if(Status >= Computed) return;
    if(Complex)
        computeImpl<true>();
    else
        computeImpl<false>();
    Status = Computed;
}

template<bool Complex> void HamiltonianPart::computeImpl()
{
    auto & HMatrix_ = *std::static_pointer_cast<MatrixType<Complex>>(HMatrix);
    if (HMatrix_.rows() == 1) {
        assert (std::abs(HMatrix_(0,0) - std::real(HMatrix_(0,0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << std::real(HMatrix_(0,0));
        HMatrix_(0,0) = 1;
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixType<Complex>> Solver(HMatrix_, Eigen::ComputeEigenvectors);
        HMatrix_ = Solver.eigenvectors();
        Eigenvalues = Solver.eigenvalues(); // eigenvectors are ready
    }
}

template<bool Complex> const MatrixType<Complex>& HamiltonianPart::getMatrix() const {
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<const MatrixType<Complex>>(HMatrix);
}
template const MatrixType<true>& HamiltonianPart::getMatrix<true>() const;
template const MatrixType<false>& HamiltonianPart::getMatrix<false>() const;

template<bool Complex> MatrixType<Complex>& HamiltonianPart::getMatrix() {
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return *std::static_pointer_cast<MatrixType<Complex>>(HMatrix);
}
template MatrixType<true>& HamiltonianPart::getMatrix<true>();
template MatrixType<false>& HamiltonianPart::getMatrix<false>();

RealType HamiltonianPart::getEigenValue(InnerQuantumState state) const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues(state);
}

const RealVectorType& HamiltonianPart::getEigenValues() const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues;
}

InnerQuantumState HamiltonianPart::getSize() const
{
    return S.getBlockSize(Block);
}

void HamiltonianPart::print_to_screen() const
{
    if(Complex)
        INFO(*std::static_pointer_cast<MatrixType<true>>(HMatrix) << std::endl);
    else
        INFO(*std::static_pointer_cast<MatrixType<false>>(HMatrix) << std::endl);
}

template<bool Complex>
VectorType<Complex> HamiltonianPart::getEigenState(InnerQuantumState state) const
{
    if (Status < Computed) throw exStatusMismatch();
    if(Complex != isComplex())
        throw std::runtime_error("Stored matrix type mismatch (real/complex)");
    return std::static_pointer_cast<MatrixType<Complex>>(HMatrix)->col(state);
}

RealType HamiltonianPart::getMinimumEigenvalue() const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues.minCoeff();
}

bool HamiltonianPart::reduce(RealType ActualCutoff)
{
    if (Status < Computed) throw exStatusMismatch();
    InnerQuantumState counter=0;
    for (counter=0; (counter< (unsigned int)Eigenvalues.size() && Eigenvalues[counter]<=ActualCutoff); ++counter){};
    std::cout << "Left " << counter << " eigenvalues : " << std::endl;

    if (counter) {
        std::cout << Eigenvalues.head(counter) << std::endl << "_________" << std::endl;
        Eigenvalues = Eigenvalues.head(counter);
        if(Complex) {
            auto & HMatrix_ = *std::static_pointer_cast<MatrixType<true>>(HMatrix);
            HMatrix_ = HMatrix_.topLeftCorner(counter,counter);
        } else {
            auto & HMatrix_ = *std::static_pointer_cast<MatrixType<false>>(HMatrix);
            HMatrix_ = HMatrix_.topLeftCorner(counter,counter);
        }
      return true;
    } else
        return false;
}

#ifdef ENABLE_SAVE_PLAINTEXT
bool HamiltonianPart<Complex>::savetxt(const boost::filesystem::path &path1)
{
    boost::filesystem::create_directory(path1);
    boost::filesystem::fstream out;
    if (Status >= Computed) {
        out.open(path1 / boost::filesystem::path("evals.dat"),std::ios_base::out);
        out << Eigenvalues << std::endl;
        out.close();
        out.open(path1 / boost::filesystem::path("evals_shift.dat"),std::ios_base::out);
        out << __num_format<RealVectorType>(Eigenvalues - RealMatrixType::Identity(Eigenvalues.size(),Eigenvalues.size()).diagonal()*getMinimumEigenvalue()) << std::endl;
        out.close();
        };
    if (Status >= Prepared) {
        out.open(path1 / boost::filesystem::path("evecs.dat"),std::ios_base::out);
        out << H << std::endl;
        out.close();
        };
    out.open(path1 / boost::filesystem::path("info.dat"),std::ios_base::out);
    out << "Quantum numbers: " << QN << std::endl;
    out << "Block number:    " << Block << std::endl;
    out.close();
    return true;
}
#endif

} // end of namespace Pomerol

