/** \file include/pomerol/Hamiltonian.h
** \brief A storage of the matrix elements of the hamiltonian in Fock basis, provides eigenvalues and eigenfunctions
**
** \author Andrey Antipov(Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_HAMILTONIAN_H
#define __INCLUDE_HAMILTONIAN_H

#include <memory>

#include <mpi.h>

#include "Misc.h"
#include "IndexClassification.h"
#include "IndexHamiltonian.h"
#include "StatesClassification.h"
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
template<bool Complex = false>
class Hamiltonian : public ComputableObject
{

public:

    using PartT = HamiltonianPart<Complex>;

private:
    /** Array of pointers to the Hamiltonian Parts */
    std::vector<std::shared_ptr<PartT> > parts;
    /** A reference to the IndexClassification object. */
    const IndexClassification<Complex> &IndexInfo;
    /** A reference to the IndexHamiltonian object. */
    const IndexHamiltonian<Complex> &F;
    /** A reference to the StatesClassification object. */
    const StatesClassification<Complex>& S;
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy;
public:

    /** Constructor. */
    Hamiltonian(const IndexClassification<Complex> &IndexInfo,
                const IndexHamiltonian<Complex>& F,
                const StatesClassification<Complex> &S);
    /** Destructor. */
    ~Hamiltonian();

    void prepare(const MPI_Comm &comm = MPI_COMM_WORLD);
    void compute(const MPI_Comm &comm = MPI_COMM_WORLD);
    void reduce(const RealType Cutoff);

    const PartT& getPart(const QuantumNumbers<Complex> &in) const;
    const PartT& getPart(BlockNumber in) const;
    RealType getEigenValue(unsigned long state) const;
    RealVectorType getEigenValues() const;
    RealType getGroundEnergy() const;

    /** Save the data to the directory.
     * \param[in] path Path to the directory.
     */
    #ifdef ENABLE_SAVE_PLAINTEXT
    bool savetxt(const boost::filesystem::path &path);
    #endif

private:
    void computeGroundEnergy();
};

extern template class Hamiltonian<false>;
extern template class Hamiltonian<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIAN_H

