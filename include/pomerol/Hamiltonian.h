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
class Hamiltonian : public ComputableObject
{
    /** Array of pointers to the Hamiltonian Parts */
    std::vector<std::shared_ptr<HamiltonianPart> > parts;
    /** A reference to the IndexClassification object. */
    const IndexClassification &IndexInfo;
    /** A reference to the IndexHamiltonian object. */
    const IndexHamiltonian &F;
    /** A reference to the StatesClassification object. */
    const StatesClassification& S;
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy;
public:

    /** Constructor. */
    Hamiltonian(const IndexClassification &IndexInfo, const IndexHamiltonian& F, const StatesClassification &S);
    /** Destructor. */
    ~Hamiltonian();

    void prepare(const MPI_Comm &comm = MPI_COMM_WORLD);
    void compute(const MPI_Comm &comm = MPI_COMM_WORLD);
    void reduce(const RealType Cutoff);

    const HamiltonianPart& getPart(const QuantumNumbers &in) const;
    const HamiltonianPart& getPart(BlockNumber in) const;
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

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIAN_H

