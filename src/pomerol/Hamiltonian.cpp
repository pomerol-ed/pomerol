#include "pomerol/Hamiltonian.h"
#include "mpi_dispatcher/misc.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#endif

namespace Pomerol{

Hamiltonian::Hamiltonian(const IndexClassification &IndexInfo, const IndexHamiltonian& F, const StatesClassification &S):
    ComputableObject(), IndexInfo(IndexInfo), F(F), S(S)
{}

Hamiltonian::~Hamiltonian()
{
}

void Hamiltonian::prepare(const MPI_Comm& comm)
{
    if (Status >= Prepared) return;
    BlockNumber NumberOfBlocks = S.NumberOfBlocks();
    parts.resize(NumberOfBlocks);
    int rank = pMPI::rank(comm);
    if (!rank) INFO_NONEWLINE("Preparing Hamiltonian parts...");


    for (BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; CurrentBlock++)
    {
	    parts[CurrentBlock].reset(new HamiltonianPart(IndexInfo,F, S, CurrentBlock));
    }
    pMPI::mpi_skel<pMPI::PrepareWrap<HamiltonianPart> > skel;
    skel.parts.resize(parts.size());
    for (size_t i=0; i<parts.size(); i++) { skel.parts[i] = pMPI::PrepareWrap<HamiltonianPart>(*parts[i]);};
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm,false);
    MPI_Barrier(comm);
    for (size_t p = 0; p<parts.size(); p++) {
        if (rank == job_map[p]){
            if (parts[p]->Status != HamiltonianPart::Prepared) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }

            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      MPI_MELEM_DATATYPE,
                      rank,
                      comm
                     );
        } else {
            auto mat_rows = parts[p]->getSize();
            parts[p]->H.resize(mat_rows, mat_rows);
            MPI_Bcast(parts[p]->H.data(),
                      mat_rows * mat_rows,
                      MPI_MELEM_DATATYPE,
                      job_map[p],
                      comm
                     );
            parts[p]->Status = HamiltonianPart::Prepared;
        }
    }
    Status = Prepared;
}

void Hamiltonian::compute(const MPI_Comm& comm)
{
    if (Status >= Computed) return;

    // Create a "skeleton" class with pointers to part that can call a compute method
    pMPI::mpi_skel<pMPI::ComputeWrap<HamiltonianPart> > skel;
    skel.parts.resize(parts.size());
    for (size_t i = 0; i < parts.size(); i++) {
      skel.parts[i] = pMPI::ComputeWrap<HamiltonianPart>(*parts[i],parts[i]->getSize());
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm, true);
    int rank = pMPI::rank(comm);
    int comm_size = pMPI::size(comm);

    // Start distributing data
    MPI_Barrier(comm);
    for (size_t p = 0; p < parts.size(); p++) {
        if (rank == job_map[p]){
            if (parts[p]->Status != HamiltonianPart::Computed) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }
            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      MPI_MELEM_DATATYPE,
                      rank,
                      comm
                    );
            MPI_Bcast(parts[p]->Eigenvalues.data(),
                      parts[p]->H.rows(),
                      MPI_DOUBLE,
                      rank,
                      comm
                    );
        } else {
            parts[p]->Eigenvalues.resize(parts[p]->H.rows());
            MPI_Bcast(parts[p]->H.data(),
                      parts[p]->H.size(),
                      MPI_MELEM_DATATYPE,
                      job_map[p],
                      comm);
            MPI_Bcast(parts[p]->Eigenvalues.data(),
                      parts[p]->H.rows(),
                      MPI_DOUBLE,
                      job_map[p],
                      comm
                     );
            parts[p]->Status = HamiltonianPart::Computed;
        }
    }

    computeGroundEnergy();
    Status = Computed;
}

void Hamiltonian::reduce(const RealType Cutoff)
{
    std::cout << "Performing EV cutoff at " << Cutoff << " level" << std::endl;
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++)
    {
	parts[CurrentBlock]->reduce(GroundEnergy+Cutoff);
    }
}

void Hamiltonian::computeGroundEnergy()
{
    RealVectorType LEV(size_t(S.NumberOfBlocks()));
    BlockNumber NumberOfBlocks = parts.size();
    for (BlockNumber CurrentBlock=0; CurrentBlock<NumberOfBlocks; CurrentBlock++) {
	    LEV(static_cast<int>(CurrentBlock),0) = parts[CurrentBlock]->getMinimumEigenvalue();
    }
    GroundEnergy=LEV.minCoeff();
}

const HamiltonianPart& Hamiltonian::getPart(const QuantumNumbers &in) const
{
    return *parts[S.getBlockNumber(in)];
}

const HamiltonianPart& Hamiltonian::getPart(BlockNumber in) const
{
    return *parts[in];
}

RealType Hamiltonian::getEigenValue(QuantumState state) const
{
    InnerQuantumState InnerState = S.getInnerState(state);
    return getPart(S.getBlockNumber(state)).getEigenValue(InnerState);
}

RealVectorType Hamiltonian::getEigenValues() const
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


RealType Hamiltonian::getGroundEnergy() const
{
    return GroundEnergy;
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
