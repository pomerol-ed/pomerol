#include "pomerol/Hamiltonian.h"
#include "mpi_dispatcher/misc.hpp"

#include <algorithm>

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#endif

namespace Pomerol {

template<bool Complex>
void Hamiltonian::compute_impl(const MPI_Comm& comm) {

    using PartType = HamiltonianPart<Complex>;

    // Create a "skeleton" class with pointers to part that can call a compute method
    pMPI::mpi_skel<pMPI::ComputeWrap<PartType>> skel;
    skel.parts.resize(parts.size());
    for (size_t p = 0; p < parts.size(); p++) {
      auto part = std::static_pointer_cast<PartType>(parts[p]);
      skel.parts[p] = pMPI::ComputeWrap<PartType>(*part, part->getSize());
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true);
    int rank = pMPI::rank(comm);

    // Start distributing data
    MPI_Barrier(comm);
    MPI_Datatype H_dt = Complex ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;
    for (size_t p = 0; p < parts.size(); ++p) {
        auto part = std::static_pointer_cast<PartType>(parts[p]);
        if (rank == job_map[p]){
            if (part->Status != PartType::Computed) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            MPI_Bcast(part->H.data(), part->H.size(), H_dt, rank, comm);
            MPI_Bcast(part->Eigenvalues.data(),
                      part->Eigenvalues.size(),
                      MPI_DOUBLE,
                      rank,
                      comm
                    );
        } else {
            part->Eigenvalues.resize(part->H.rows());
            MPI_Bcast(part->H.data(), part->H.size(), H_dt, job_map[p], comm);
            MPI_Bcast(part->Eigenvalues.data(),
                      part->Eigenvalues.size(),
                      MPI_DOUBLE,
                      job_map[p],
                      comm
                     );
            part->Status = PartType::Computed;
        }
    }
}

void Hamiltonian::compute(const MPI_Comm& comm)
{
    if (Status >= Computed) return;

    if(Complex)
        compute_impl<true>(comm);
    else
        compute_impl<false>(comm);

    computeGroundEnergy();
    Status = Computed;
}

void Hamiltonian::reduce(const RealType Cutoff)
{
    std::cout << "Performing EV cutoff at " << Cutoff << " level" << std::endl;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; ++CurrentBlock)
    {
        if(Complex)
            getComplexPart(CurrentBlock).reduce(GroundEnergy+Cutoff);
        else
            getRealPart(CurrentBlock).reduce(GroundEnergy+Cutoff);
    }
}

void Hamiltonian::computeGroundEnergy()
{
    RealVectorType LEV(size_t(S.getNumberOfBlocks()));
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; ++CurrentBlock) {
        LEV(static_cast<int>(CurrentBlock), 0) = Complex ?
            getComplexPart(CurrentBlock).getMinimumEigenvalue() :
            getRealPart(CurrentBlock).getMinimumEigenvalue();
    }
    GroundEnergy = LEV.minCoeff();
}

RealType Hamiltonian::getEigenValue(QuantumState state) const
{
    InnerQuantumState InnerState = S.getInnerState(state);
    BlockNumber Block = S.getBlockNumber(state);
    return Complex ? getComplexPart(Block).getEigenValue(InnerState) :
                     getRealPart(Block).getEigenValue(InnerState);
}

RealVectorType Hamiltonian::getEigenValues() const
{
    RealVectorType out(S.getNumberOfStates());
    size_t copied_size = 0;
    for (BlockNumber CurrentBlock=0; CurrentBlock < S.getNumberOfBlocks(); CurrentBlock++) {
        const RealVectorType& tmp = Complex ?
            getComplexPart(CurrentBlock).getEigenValues() :
            getRealPart(CurrentBlock).getEigenValues();
        std::copy(tmp.data(), tmp.data() + tmp.size(), out.data() + copied_size);
        copied_size += tmp.size();
    }
    return out;
}

#ifdef ENABLE_SAVE_PLAINTEXT
bool Hamiltonian::savetxt(const boost::filesystem::path &path)
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

} // end of namespace Pomerol
