#include "pomerol/LibcommuteEigen.h"
#include "pomerol/HamiltonianPart.h"

#include <libcommute/loperator/mapped_basis_view.hpp>

#include <Eigen/Eigenvalues>

#include <limits>
#include <sstream>

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#endif

//
// class HamiltonianPart
//

namespace Pomerol {

template<bool Complex>
void HamiltonianPart<Complex>::prepare()
{
    if(Status >= Prepared) return;

    auto BlockSize = S.getBlockSize(Block);
    H.resize(BlockSize, BlockSize);

    auto mapper = libcommute::basis_mapper(S.getFockStates(Block));

    VectorType<Complex> ket = VectorType<Complex>::Zero(BlockSize);
    auto ket_view = mapper.make_const_view(ket);

    for(InnerQuantumState st = 0; st < BlockSize; ++st) {
        auto bra_view = mapper.make_view_no_ref(H.col(st));
        ket(st) = 1.0;
        HOp(ket_view, bra_view);
        ket(st) = .0;
    }

    assert((H.adjoint() - H).array().abs().maxCoeff()
           < 100*std::numeric_limits<RealType>::epsilon());

    Status = Prepared;
}

template<bool Complex>
void HamiltonianPart<Complex>::compute()		//method of diagonalization classificated part of Hamiltonian
{
    if (Status >= Computed) return;

    if (H.rows() == 1) {
        assert (std::abs(H(0,0) - std::real(H(0,0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << std::real(H(0,0));
        H(0,0) = 1;
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixType<Complex>> Solver(H,Eigen::ComputeEigenvectors);
        H = Solver.eigenvectors();
        Eigenvalues = Solver.eigenvalues();	// eigenvectors are ready
    }
    Status = Computed;
}

template<bool Complex>
MelemType<Complex> HamiltonianPart<Complex>::getMatrixElement(InnerQuantumState m, InnerQuantumState n) const
{
    return H(m,n);
}

template<bool Complex>
RealType HamiltonianPart<Complex>::getEigenValue(InnerQuantumState state) const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues(state);
}

template<bool Complex>
const RealVectorType& HamiltonianPart<Complex>::getEigenValues() const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues;
}

template<bool Complex>
InnerQuantumState HamiltonianPart<Complex>::getSize(void) const
{
    return S.getBlockSize(Block);
}

template<bool Complex>
void HamiltonianPart<Complex>::print_to_screen() const
{
    INFO(H << std::endl);
}

template<bool Complex>
VectorType<Complex> HamiltonianPart<Complex>::getEigenState(InnerQuantumState state) const
{
    if (Status < Computed) throw exStatusMismatch();
    return H.col(state);
}

template<bool Complex>
RealType HamiltonianPart<Complex>::getMinimumEigenvalue() const
{
    if (Status < Computed) throw exStatusMismatch();
    return Eigenvalues.minCoeff();
}

template<bool Complex>
bool HamiltonianPart<Complex>::reduce(RealType ActualCutoff)
{
    if (Status < Computed) throw exStatusMismatch();
    InnerQuantumState counter=0;
    for (counter=0; (counter< (unsigned int)Eigenvalues.size() && Eigenvalues[counter]<=ActualCutoff); ++counter){};
    std::cout << "Left " << counter << " eigenvalues : " << std::endl;
    if (counter)
	{std::cout << Eigenvalues.head(counter) << std::endl << "_________" << std::endl;
	Eigenvalues = Eigenvalues.head(counter);
	H = H.topLeftCorner(counter,counter);
	return true;
    }
    else return false;
}

#ifdef ENABLE_SAVE_PLAINTEXT
template<bool Complex>
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

template class HamiltonianPart<false>;
template class HamiltonianPart<true>;

} // end of namespace Pomerol

