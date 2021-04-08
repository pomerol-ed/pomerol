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
    struct {
        RealType ActualCutoff;
        void operator()(BlockNumber, RealPartType & p) { p.reduce(ActualCutoff); }
        void operator()(BlockNumber, ComplexPartType & p) { p.reduce(ActualCutoff); }
    } visitor;
    visitor.ActualCutoff = GroundEnergy + Cutoff;
    visitAllParts(visitor);
}

InnerQuantumState Hamiltonian::getBlockSize(BlockNumber Block) const
{
    struct {
        InnerQuantumState operator()(RealPartType const& p) { return p.getSize(); }
        InnerQuantumState operator()(ComplexPartType const& p) { return p.getSize(); }
    } visitor;
    return visitPart(Block, visitor);
}

void Hamiltonian::computeGroundEnergy()
{
    struct {
        RealVectorType lev;
        void operator()(BlockNumber b, RealPartType & p) { lev(b) = p.getMinimumEigenvalue(); }
        void operator()(BlockNumber b, ComplexPartType & p) { lev(b) = p.getMinimumEigenvalue(); }
    } visitor;
    visitor.lev.resize(S.getNumberOfBlocks());
    visitAllParts(visitor);
    GroundEnergy = visitor.lev.minCoeff();
}

RealType Hamiltonian::getEigenValue(QuantumState state) const
{
    struct {
        InnerQuantumState InnerState;
        RealType operator()(RealPartType const& p) { return p.getEigenValue(InnerState); }
        RealType operator()(ComplexPartType const& p) { return p.getEigenValue(InnerState); }
    } visitor {S.getInnerState(state)};
    return visitPart(S.getBlockNumber(state), visitor);
}

RealVectorType const& Hamiltonian::getEigenValues(BlockNumber Block) const
{
    struct {
        RealVectorType const& operator()(RealPartType const& p) { return p.getEigenValues(); }
        RealVectorType const& operator()(ComplexPartType const& p) { return p.getEigenValues(); }
    } visitor;
    return visitPart(Block, visitor);
}

RealVectorType Hamiltonian::getEigenValues() const
{
    struct {
        RealVectorType out;
        size_t copied_size = 0;
        void add_ev_chunk(RealVectorType const& ev) {
            out.segment(copied_size, ev.size()) = ev;
            copied_size += ev.size();
        }
        void operator()(BlockNumber, RealPartType const& p) { add_ev_chunk(p.getEigenValues()); }
        void operator()(BlockNumber, ComplexPartType const& p) { add_ev_chunk(p.getEigenValues()); }
    } visitor;
    visitor.out.resize(S.getNumberOfStates());
    visitAllParts(visitor);
    return visitor.out;
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
