//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


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
    /** A value of the ground energy - needed for further renormalization */
    RealType GroundEnergy;
public:

    Hamiltonian(IndexClassification &F_, StatesClassification &S_);
    ~Hamiltonian();

    void prepare(); // was void enter();

    HamiltonianPart& part(const QuantumNumbers &in);
    HamiltonianPart& part(BlockNumber in);
    RealType eigenval( QuantumState &state );
    RealType getGroundEnergy();

    void compute(); // was void diagonalize();

    void dump();
    void reduce(const RealType Cutoff);

    void save(H5::CommonFG* RootGroup) const;
    void load(const H5::CommonFG* RootGroup);

private:
    void computeGroundEnergy();
};

#endif // endif :: #ifndef __INCLUDE_HAMILTONIAN_H

