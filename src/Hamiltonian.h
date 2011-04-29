/** \file src/DensityMatrix.h
** \brief A storage of the matrix elements of the hamiltonian in Fock basis, provides eigenvalues and eigenfunctions
** 
** \author Andrey Antipov(antipov@ct-qmc.org)
** \author Igor Krivenko (igor@shg.ru)
*/

#ifndef __INCLUDE_HAMILTONIAN_H
#define __INCLUDE_HAMILTONIAN_H
#include "Misc.h"
#include "HDF5Storage.h"
#include "ComputableObject.h"
#include "IndexClassification.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "output.h"

/** This class represents a Hamiltonian, written as a matrix of matrix elements in a Fock basis.
 * It is a container for several hamiltonian parts, each for single defined QuantumNumbers and a corresponding BlockNumber. 
 * It provides eigenvalues and eigenfunctions of any of its parts once they are obtained within its parts. 
 * The diagonalization and entering routines are done inside Hamiltonian Parts
 */
class Hamiltonian : public ComputableObject, public HDF5Storable
{
    /** Array of pointers to the Hamiltonian Parts */
    std::vector<HamiltonianPart*> parts;
    /** A reference to the object, which contains all info about how sites and spins of the lattice are defined as bits */
    IndexClassification &Formula;
    /** Reference to a states classification object. */
    StatesClassification& S;
    /** Reference to an output handling object */
    output_handle &OUT;
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy;
public:

    Hamiltonian(IndexClassification &F_, StatesClassification &S_,output_handle &OUT_, std::string &config_path_);
    ~Hamiltonian();
    void enter();

    HamiltonianPart& part(const QuantumNumbers &in);
    HamiltonianPart& part(BlockNumber in);
    RealType eigenval( QuantumState &state );
    RealType getGroundEnergy();

    void diagonalize();
    void dump();
    void reduce(const RealType Cutoff);

    void save(H5::CommonFG* FG) const;
    void load(const H5::CommonFG* FG);

private:
    void computeGroundEnergy();
};

#endif // endif :: #ifndef __INCLUDE_HAMILTONIAN_H

