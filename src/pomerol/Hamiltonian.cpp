#include "pomerol/Hamiltonian.h"
#include "mpi_dispatcher/misc.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#endif

namespace Pomerol{

template<bool Complex>
Hamiltonian<Complex>::Hamiltonian(const IndexClassification<Complex> &IndexInfo,
                                  const IndexHamiltonian<Complex>& F,
                                  const StatesClassification<Complex> &S):
    ComputableObject(), IndexInfo(IndexInfo), F(F), S(S)
{}

template<bool Complex>
Hamiltonian<Complex>::~Hamiltonian()
{
}

template<bool Complex>
void Hamiltonian<Complex>::prepare(const MPI_Comm& comm)
{
    if (Status >= Prepared) return;
    BlockNumber NumberOfBlocks = S.NumberOfBlocks();
    parts.resize(NumberOfBlocks);
    int rank = pMPI::rank(comm);
    if (!rank) INFO_NONEWLINE("Preparing Hamiltonian parts...");


    for (BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; CurrentBlock++)
    {
	    parts[CurrentBlock].reset(new PartT(IndexInfo,F, S, CurrentBlock));
    }
    pMPI::mpi_skel<pMPI::PrepareWrap<PartT> > skel;
    skel.parts.resize(parts.size());
    for (size_t i=0; i<parts.size(); i++) {
      skel.parts[i] = pMPI::PrepareWrap<PartT>(*parts[i]);
    };
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm,false);
    MPI_Barrier(comm);
    for (size_t p = 0; p<parts.size(); p++) {
        if (rank == job_map[p]){
            if (parts[p]->Status != PartT::Prepared) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }

            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      pMPI::mpi_datatype<MelemType<Complex>>(),
                      rank,
                      comm
                     );
        } else {
            auto mat_rows = parts[p]->getSize();
            parts[p]->H.resize(mat_rows, mat_rows);
            MPI_Bcast(parts[p]->H.data(),
                      mat_rows * mat_rows,
                      pMPI::mpi_datatype<MelemType<Complex>>(),
                      job_map[p],
                      comm
                     );
            parts[p]->Status = PartT::Prepared;
        }
    }
    Status = Prepared;
}

template<bool Complex>
void Hamiltonian<Complex>::compute(const MPI_Comm& comm)
{
    if (Status >= Computed) return;

    // Create a "skeleton" class with pointers to part that can call a compute method
    pMPI::mpi_skel<pMPI::ComputeWrap<PartT>> skel;
    skel.parts.resize(parts.size());
    for (size_t i = 0; i < parts.size(); i++) {
      skel.parts[i] = pMPI::ComputeWrap<PartT>(*parts[i],parts[i]->getSize());
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true);
    int rank = pMPI::rank(comm);

    // Start distributing data
    MPI_Barrier(comm);
    for (size_t p = 0; p < parts.size(); p++) {
        if (rank == job_map[p]){
            if (parts[p]->Status != PartT::Computed) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      pMPI::mpi_datatype<MelemType<Complex>>(),
                      rank,
                      comm
                    );
            MPI_Bcast(parts[p]->Eigenvalues.data(),
                      parts[p]->Eigenvalues.size(),
                      MPI_DOUBLE,
                      rank,
                      comm
                    );
        } else {
            parts[p]->Eigenvalues.resize(parts[p]->H.rows());
            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      pMPI::mpi_datatype<MelemType<Complex>>(),
                      job_map[p],
                      comm);
            MPI_Bcast(parts[p]->Eigenvalues.data(),
                      parts[p]->Eigenvalues.size(),
                      MPI_DOUBLE,
                      job_map[p],
                      comm
                     );
            parts[p]->Status = PartT::Computed;
        }
    }

    computeGroundEnergy();
    Status = Computed;
}

template<bool Complex>
void Hamiltonian<Complex>::reduce(const RealType Cutoff)
{
    std::cout << "Performing EV cutoff at " << Cutoff << " level" << std::endl;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++)
    {
	parts[CurrentBlock]->reduce(GroundEnergy+Cutoff);
    }
}

template<bool Complex>
void Hamiltonian<Complex>::computeGroundEnergy()
{
    RealVectorType LEV(size_t(S.NumberOfBlocks()));
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++) {
	    LEV(static_cast<int>(CurrentBlock),0) = parts[CurrentBlock]->getMinimumEigenvalue();
    }
    GroundEnergy=LEV.minCoeff();
}

template<bool Complex>
auto Hamiltonian<Complex>::getPart(const QuantumNumbers<Complex> &in) const -> const PartT&
{
    return *parts[S.getBlockNumber(in)];
}

template<bool Complex>
auto Hamiltonian<Complex>::getPart(BlockNumber in) const -> const PartT&
{
    return *parts[in];
}

template<bool Complex>
RealType Hamiltonian<Complex>::getEigenValue(QuantumState state) const
{
    InnerQuantumState InnerState = S.getInnerState(state);
    return getPart(S.getBlockNumber(state)).getEigenValue(InnerState);
}

template<bool Complex>
RealVectorType Hamiltonian<Complex>::getEigenValues() const
{
    RealVectorType out(S.getNumberOfStates());
    size_t i=0;
    for (BlockNumber CurrentBlock=0; CurrentBlock<S.NumberOfBlocks(); CurrentBlock++) {
        const RealVectorType& tmp = parts[CurrentBlock]->getEigenValues();
        std::copy(tmp.data(), tmp.data() + tmp.size(), out.data()+i);
        i+=tmp.size();
    }
    return out;
}

template<bool Complex>
RealType Hamiltonian<Complex>::getGroundEnergy() const
{
    return GroundEnergy;
}

#ifdef ENABLE_SAVE_PLAINTEXT
template<bool Complex>
bool Hamiltonian<Complex>::savetxt(const boost::filesystem::path &path)
{
    BlockNumber NumberOfBlocks = parts.size();
    boost::filesystem::create_directory(path);
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++) {
        std::stringstream tmp;
        tmp << "part" << S.getQuantumNumbers(CurrentBlock);
        boost::filesystem::path out = path / boost::filesystem::path (tmp.str());
	    parts[CurrentBlock]->savetxt(out);
        }
    return true;
}
#endif

template class Hamiltonian<false>;
template class Hamiltonian<true>;

} // end of namespace Pomerol
