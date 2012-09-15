//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
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

/** \file src/HamiltonianPart.h
** \brief Declaration of the HamiltonianPart class.
**
** \author Andrey Antipov (antipov@ct-qmc.org)
** \author Igor Krivenko (igor@shg.ru)
*/

#ifndef __INCLUDE_HAMILTONIANPART_H
#define __INCLUDE_HAMILTONIANPART_H

#include "Misc.h"
#include "HDF5Storage.h"
#include "IndexClassification.h"
#include "StatesClassification.h"

namespace Pomerol{

/** HamiltonianPart is a class, which stores and diagonalizes the block of the Hamiltonian, which corresponds to a set of given quantum numbers. */
class HamiltonianPart : public HDF5Storable, public ComputableObject {
    /** Computation statuses of the object. */
    enum {Constructed, Prepared, Diagonalized};

    /** A reference to the IndexClassification object. */
    const IndexClassification &IndexInfo;
    /** A reference to the IndexHamiltonian object. */
    const IndexHamiltonian &F;
    /** A reference to the StateClassification object. */
    const StatesClassification &S;

    /** The number of Block. Defined in StatesClassification. */
    BlockNumber Block;
    /** QuantumNumbers of the block. Consructed in Symmetrizer and defined in StatesClassification. */
    QuantumNumbers QN;

    /** A matrix filled with matrix elements of HamiltonianPart in the space of FockState's.
     *  After diagonalization it stores the eigenfunctions of the problem in a rows of H. */
    MatrixType H;                
    /** A vector of eigenvalues of the HamiltonianPart. */
    RealVectorType Eigenvalues;      

public:

    /** Constructor.
     * \param[in] IndexInfo IndexClassification object. Provides information about the indices in the problem. 
     * \param[in] F IndexHamiltonian object. Provides all Terms required to fill the Hamiltonian part.
     * \param[in] S StatesClassification object. Provides information about Fock States of the problem.
     * \param[in] Block The BlockNumber of current part. It is a genuine id of the part. */
    HamiltonianPart(const IndexClassification &IndexInfo, const IndexHamiltonian &F, const StatesClassification &S, const BlockNumber& Block);

    /** Fill in the H matrix. */
    void prepare(void);
    /** Diagonalized the H matrix and get EigenValues. */
    void diagonalize(void);
    
    bool reduce(RealType ActualCutoff); // Useless now

    /** Return the total dimensionality of the H matrix. This corresponds to the one in StatesClassfication. */
    InnerQuantumState getSize(void) const;

    /** Get the matrix element of the HamiltonianPart by the number of states inside the part. */ 
    MelemType getMatrixElement(InnerQuantumState m, InnerQuantumState n) const; //return H(m,n)
    /** Get the matrix element of the Hamiltonian within two given FockStates. */
    MelemType getMatrixElement(FockState m, FockState n) const; //return H(m,n)

    /** Get the eigenvalue of the H matrix.
     * \param[in] Number of eigenvalue. */
    RealType getEigenValue(InnerQuantumState state) const; 

    /** Returns calculated eigenvalues. */
    const RealVectorType& getEigenValues() const; 

    /** Return the hamiltonian part matrix. */
    const MatrixType& getMatrix() const;

    /** Return the lowest Eigenvalue of the current part. */
    RealType getMinimumEigenvalue() const;        
    /** Return the eigenstate of the H matrix.
     * \param[in] Number of eigenvalue. */
    VectorType getEigenState(InnerQuantumState state) const;

    /** Return the QuantumNumbers associated with the Hamiltonian part. */
    QuantumNumbers getQuantumNumbers() const; 
    /** Return the BlockNumber associated with the Hamiltonian part. */
    BlockNumber getBlockNumber() const;

    /** Print the part of hamiltonian to screen. */
    void print_to_screen() const; 

    #ifdef POMEROL_USE_PLAIN_SAVE
    /** Save the data to the file.
     * \param[in] path Path to the file.
     */
    bool savetxt(const boost::filesystem::path &path);
    #endif
    
    /** Save the HamiltonianPart to the HDF5 group. */
    void save(H5::CommonFG* RootGroup) const;
    /** Load the HamiltonianPart from the HDF5 group. */
    void load(const H5::CommonFG* RootGroup);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIANPART_H
