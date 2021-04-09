/** \file include/pomerol/Hamiltonian.h
** \brief A storage of the matrix elements of the hamiltonian in Fock basis, provides eigenvalues and eigenfunctions
**
** \author Andrey Antipov(Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_HAMILTONIAN_H
#define __INCLUDE_HAMILTONIAN_H

#include <memory>
#include <type_traits>

#include "mpi_dispatcher/misc.hpp"
#include "mpi_dispatcher/mpi_skel.hpp"

#include "Misc.h"
#include "Operators.h"
#include "IndexClassification.h"
#include "StatesClassification.h"
#include "HilbertSpace.h"
#include "HamiltonianPart.h"

#ifdef ENABLE_SAVE_PLAINTEXT
#include <boost/filesystem/path.hpp>
#endif

namespace Pomerol{

/** This class represents a Hamiltonian, written as a matrix of matrix elements in a Fock basis.
 * It is a container for several hamiltonian parts, each for single defined QuantumNumbers and a corresponding BlockNumber.
 * It provides eigenvalues and eigenfunctions of any of its parts once they are obtained within its parts.
 * The diagonalization and entering routines are done inside of HamiltonianPart instances.
 */
class Hamiltonian : public ComputableObject
{
    bool Complex;

    /** Array of pointers to the Hamiltonian Parts */
    std::vector<std::shared_ptr<void>> parts;
    /** A reference to the StatesClassification object. */
    const StatesClassification& S;
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy;
public:

    using RealPartType = HamiltonianPart<false>;
    using ComplexPartType = HamiltonianPart<true>;

    /** Constructor. */
    Hamiltonian(const StatesClassification &S) : S(S) {}

    template<typename ScalarType, typename... IndexTypes>
    void prepare(const Operators::expression<ScalarType, IndexTypes...> &H,
                 const HilbertSpace<ScalarType, IndexTypes...> &Symm,
                 const MPI_Comm &comm = MPI_COMM_WORLD);
    void compute(const MPI_Comm &comm = MPI_COMM_WORLD);
    void reduce(const RealType Cutoff);

    bool isComplex() const { return Complex; }

    InnerQuantumState getBlockSize(BlockNumber Block) const;

    template<typename F>
    auto visitPart(BlockNumber Block, F && f) ->
        typename std::result_of<F(RealPartType &)>::type
    {
        if(Complex)
            return f(*std::static_pointer_cast<ComplexPartType>(parts[Block]));
        else
            return f(*std::static_pointer_cast<RealPartType>(parts[Block]));
    }
    template<typename F>
    auto visitPart(BlockNumber Block, F && f) const ->
        typename std::result_of<F(RealPartType const&)>::type
    {
        if(Complex)
            return f(*std::static_pointer_cast<const ComplexPartType>(parts[Block]));
        else
            return f(*std::static_pointer_cast<const RealPartType>(parts[Block]));
    }

    template<typename F>
    void visitAllParts(F &&f)
    {
        for(BlockNumber Block = 0; Block < parts.size(); ++Block) {
            if(Complex)
                f(Block, *std::static_pointer_cast<ComplexPartType>(parts[Block]));
            else
                f(Block, *std::static_pointer_cast<RealPartType>(parts[Block]));
        }
    }
    template<typename F>
    void visitAllParts(F &&f) const
    {
        for(BlockNumber Block = 0; Block < parts.size(); ++Block) {
            if(Complex)
                f(Block, *std::static_pointer_cast<const ComplexPartType>(parts[Block]));
            else
                f(Block, *std::static_pointer_cast<const RealPartType>(parts[Block]));
        }
    }

    RealType getEigenValue(unsigned long state) const;
    RealVectorType const& getEigenValues(BlockNumber Block) const;
    RealVectorType getEigenValues() const;
    RealType getGroundEnergy() const { return GroundEnergy; }

    /** Save the data to the directory.
     * \param[in] path Path to the directory.
     */
    #ifdef ENABLE_SAVE_PLAINTEXT
    bool savetxt(const boost::filesystem::path &path);
    #endif

private:
    void computeGroundEnergy();

    template<bool C> void compute_impl(const MPI_Comm& comm);
};

template<typename ScalarType,
         typename... IndexTypes>
void Hamiltonian::prepare(Operators::expression<ScalarType, IndexTypes...> const& H,
                          const HilbertSpace<ScalarType, IndexTypes...> &HS,
                          const MPI_Comm &comm) {

    if (Status >= Prepared) return;

    constexpr bool C = std::is_same<ScalarType, ComplexType>::value;
    Complex = C;

    BlockNumber NumberOfBlocks = S.getNumberOfBlocks();
    parts.resize(NumberOfBlocks);
    int rank = pMPI::rank(comm);
    if (!rank) INFO_NONEWLINE("Preparing Hamiltonian parts...");

    auto const& FullHilbertSpace = HS.getFullHilbertSpace();

    libcommute::loperator<ScalarType, libcommute::fermion> HOp(H, FullHilbertSpace);
    for (BlockNumber CurrentBlock = 0; CurrentBlock < NumberOfBlocks; ++CurrentBlock) {
        parts[CurrentBlock].reset(new HamiltonianPart<C>(HOp, S, CurrentBlock));
    }

    using PartType = HamiltonianPart<C>;

    pMPI::mpi_skel<pMPI::PrepareWrap<PartType> > skel;
    skel.parts.resize(parts.size());
    for (size_t p = 0; p < parts.size(); ++p) {
        skel.parts[p] = pMPI::PrepareWrap<PartType>(
            *std::static_pointer_cast<PartType>(parts[p])
        );
    }
    std::map<pMPI::JobId, pMPI::WorkerId> job_map = skel.run(comm,false);
    MPI_Barrier(comm);

    MPI_Datatype H_dt = Complex ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;

    for (size_t p = 0; p < parts.size(); ++p) {
        auto part = std::static_pointer_cast<PartType>(parts[p]);

        if (rank == job_map[p]) {

            if (part->Status != HamiltonianPart<C>::Prepared) {
                ERROR ("Worker" << rank << " didn't calculate part" << p);
                throw std::logic_error("Worker didn't calculate this part.");
            }

            MPI_Bcast(part->H.data(), part->H.size(), H_dt, rank, comm);
        } else {
            auto mat_rows = part->getSize();
            part->H.resize(mat_rows, mat_rows);
            MPI_Bcast(part->H.data(), mat_rows * mat_rows, H_dt, job_map[p], comm);
            part->Status = HamiltonianPart<C>::Prepared;
        }
    }
    Status = Prepared;
}

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIAN_H

