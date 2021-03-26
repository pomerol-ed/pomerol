#include"pomerol/HamiltonianPart.h"
#include"pomerol/StatesClassification.h"
#include<sstream>
#include<Eigen/Eigenvalues>

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#endif

// class HamiltonianPart

namespace Pomerol{

template<bool Complex>
HamiltonianPart<Complex>::HamiltonianPart(const IndexClassification<Complex>& IndexInfo,
                                          const IndexHamiltonian<Complex> &F,
                                          const StatesClassification<Complex> &S,
                                          const BlockNumber& Block):
    ComputableObject(),
    IndexInfo(IndexInfo),
    F(F), S(S),
    Block(Block), QN(S.getQuantumNumbers(Block))
{
}

template<bool Complex>
void HamiltonianPart<Complex>::prepare()
{
    size_t BlockSize = S.getBlockSize(Block);

    H.resize(BlockSize,BlockSize);
    H.setZero();
    typename std::map<FockState, MelemT>::const_iterator melem_it;

    for(InnerQuantumState right_st=0; right_st<BlockSize; right_st++)
    {
        FockState ket = S.getFockState(Block,right_st);
        std::map<FockState, MelemT> mapStates = F.actRight(ket);
        for (melem_it=mapStates.begin(); melem_it!=mapStates.end(); melem_it++) {
            FockState bra = melem_it -> first;
            MelemT melem = melem_it -> second;
            //DEBUG("<" << bra << "|" << melem << "|" << F << "|" << ket << ">");
            InnerQuantumState left_st = S.getInnerState(bra);
//            if (left_st > right_st) { ERROR("!"); exit(1); };
            H(left_st,right_st) = melem;
        }
    }

//    H.triangularView<Eigen::Lower>() = H.triangularView<Eigen::Upper>().transpose();
//    assert(MatrixType(H.triangularView<Eigen::Lower>()) == MatrixType(H.triangularView<Eigen::Upper>().transpose()));
    assert((H.adjoint() - H).array().abs().maxCoeff() < 100*std::numeric_limits<RealType>::epsilon());
    Status = Prepared;
}

template<bool Complex>
void HamiltonianPart<Complex>::compute()		//method of diagonalization classificated part of Hamiltonian
{
    if (Status >= Computed) return;
    if (H.rows() == 1) {
        assert (std::abs(H(0,0) - std::real(H(0,0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << real(H(0,0));
        H(0,0) = 1;
        }
    else {
	    Eigen::SelfAdjointEigenSolver<MatrixT> Solver(H,Eigen::ComputeEigenvectors);
	    H = Solver.eigenvectors();
	    Eigenvalues = Solver.eigenvalues();	// eigenvectors are ready
    }
    Status = Computed;
}

template<bool Complex>
auto HamiltonianPart<Complex>::getMatrixElement(InnerQuantumState m, InnerQuantumState n) const -> MelemT //return  H(m,n)
{
    return H(m,n);
}

template<bool Complex>
RealType HamiltonianPart<Complex>::getEigenValue(InnerQuantumState state) const // return Eigenvalues(state)
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues(state);
}

template<bool Complex>
const RealVectorType& HamiltonianPart<Complex>::getEigenValues() const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues;
}

template<bool Complex>
InnerQuantumState HamiltonianPart<Complex>::getSize(void) const
{
    return S.getBlockSize(Block);
}

template<bool Complex>
BlockNumber HamiltonianPart<Complex>::getBlockNumber(void) const
{
    return S.getBlockNumber(QN);
}

template<bool Complex>
QuantumNumbers<Complex> HamiltonianPart<Complex>::getQuantumNumbers(void) const
{
    return QN;
}

template<bool Complex>
void HamiltonianPart<Complex>::print_to_screen() const
{
    INFO(H << std::endl);
}

template<bool Complex>
auto HamiltonianPart<Complex>::getMatrix() const -> const MatrixT&
{
    return H;
}

template<bool Complex>
VectorType<Complex> HamiltonianPart<Complex>::getEigenState(InnerQuantumState state) const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return H.col(state);
}

template<bool Complex>
RealType HamiltonianPart<Complex>::getMinimumEigenvalue() const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues.minCoeff();
}

template<bool Complex>
bool HamiltonianPart<Complex>::reduce(RealType ActualCutoff)
{
    if ( Status < Computed ) throw (exStatusMismatch());
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

