#include "pomerol/Hamiltonian.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>

namespace Pomerol {

template<bool C>
void Hamiltonian::prepareImpl(const LOperatorTypeRC<C> &HOp, const MPI_Comm& comm)
{
    BlockNumber NumberOfBlocks = S.getNumberOfBlocks();
    int rank = pMPI::rank(comm);
    if(!rank) INFO_NONEWLINE("Preparing Hamiltonian parts...");

    parts.reserve(NumberOfBlocks);
    for(BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; ++CurrentBlock) {
        parts.emplace_back(HOp, S, CurrentBlock);
    }

    pMPI::mpi_skel<pMPI::PrepareWrap<HamiltonianPart>> skel;
    skel.parts.resize(parts.size());
    for(size_t p = 0; p < parts.size(); ++p) {
        skel.parts[p] = pMPI::PrepareWrap<HamiltonianPart>(parts[p]);
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, false);
    MPI_Barrier(comm);

    MPI_Datatype H_dt = C ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;

    for(size_t p = 0; p < parts.size(); ++p) {
        auto & part = parts[p];
        if (rank == job_map[p]) {
            if (part.getStatus() != HamiltonianPart::Prepared) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            auto & H = part.getMatrix<C>();
            MPI_Bcast(H.data(), H.size(), H_dt, rank, comm);
        } else {
            part.initHMatrix<C>();
            auto & H = part.getMatrix<C>();
            MPI_Bcast(H.data(), H.rows() * H.cols(), H_dt, job_map[p], comm);
            part.setStatus(HamiltonianPart::Prepared);
        }
    }
}

template void Hamiltonian::prepareImpl<true>(const LOperatorTypeRC<true> &p, const MPI_Comm&);
template void Hamiltonian::prepareImpl<false>(const LOperatorTypeRC<false> &, const MPI_Comm&);

template<bool C>
void Hamiltonian::computeImpl(const MPI_Comm& comm)
{
    // Create a "skeleton" class with pointers to part that can call a compute method
    pMPI::mpi_skel<pMPI::ComputeWrap<HamiltonianPart>> skel;
    skel.parts.resize(parts.size());
    for (size_t p = 0; p < parts.size(); p++) {
      auto & part = parts[p];
      skel.parts[p] = pMPI::ComputeWrap<HamiltonianPart>(part, part.getSize());
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true);
    int rank = pMPI::rank(comm);

    // Start distributing data
    MPI_Barrier(comm);
    MPI_Datatype H_dt = C ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;
    for (size_t p = 0; p < parts.size(); ++p) {
        auto & part = parts[p];
        auto & H = part.getMatrix<C>();
        if (rank == job_map[p]){
            if (part.Status != HamiltonianPart::Computed) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            MPI_Bcast(H.data(), H.size(), H_dt, rank, comm);
            MPI_Bcast(part.Eigenvalues.data(),
                      part.Eigenvalues.size(),
                      MPI_DOUBLE,
                      rank,
                      comm
                    );
        } else {
            part.Eigenvalues.resize(H.rows());
            MPI_Bcast(H.data(), H.size(), H_dt, job_map[p], comm);
            MPI_Bcast(part.Eigenvalues.data(),
                      part.Eigenvalues.size(),
                      MPI_DOUBLE,
                      job_map[p],
                      comm
                     );
            part.Status = HamiltonianPart::Computed;
        }
    }
}

void Hamiltonian::compute(const MPI_Comm& comm)
{
    if (getStatus() >= Computed) return;

    if(Complex)
        computeImpl<true>(comm);
    else
        computeImpl<false>(comm);

    computeGroundEnergy();

    setStatus(Computed);
}

void Hamiltonian::reduce(const RealType Cutoff)
{
    INFO("Performing EV cutoff at " << Cutoff << " level");
    for(auto & part : parts)
        part.reduce(GroundEnergy + Cutoff);
}

InnerQuantumState Hamiltonian::getBlockSize(BlockNumber Block) const
{
    return parts[Block].getSize();
}

void Hamiltonian::computeGroundEnergy()
{
    RealVectorType lev(S.getNumberOfBlocks());
    for(BlockNumber b = 0; b < parts.size(); ++b)
        lev(b) = parts[b].getMinimumEigenvalue();
    GroundEnergy = lev.minCoeff();
}

RealType Hamiltonian::getEigenValue(QuantumState state) const
{
    return parts[S.getBlockNumber(state)].getEigenValue(S.getInnerState(state));
}

RealVectorType const& Hamiltonian::getEigenValues(BlockNumber Block) const
{
    return parts[Block].getEigenValues();
}

RealVectorType Hamiltonian::getEigenValues() const
{
    RealVectorType out(S.getNumberOfStates());
    size_t copied_size = 0;
    for(auto const& part : parts) {
        auto const& ev = part.getEigenValues();
        out.segment(copied_size, ev.size()) = ev;
        copied_size += ev.size();
    }
    return out;
}

} // namespace Pomerol
