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


#ifndef __INCLUDE_HAMILTONIANPART_H
#define __INCLUDE_HAMILTONIANPART_H
#include "Misc.h"
#include "HDF5Storage.h"
#include "IndexClassification.h"
#include "StatesClassification.h"

#warning This class requires a major cleanup and reorganization according to the TODO list.

namespace Pomerol{

class HamiltonianPart : public HDF5Storable {

    IndexClassification &IndexInfo;
    StatesClassification &S;

    QuantumNumbers QN;

    RealMatrixType H;                //part of Hamiltonian
    RealVectorType Eigenvalues;      //vector of Eigen Values

    // DEPRECATED
    void add_nTerm(InnerQuantumState st, nTerm *N);
    // DEPRECATED
    void add_nnTerm(InnerQuantumState st, nnTerm *T);
    // DEPRECATED
    void add_spinflipTerm(InnerQuantumState st, spinflipTerm *T);
    // DEPRECATED
    int measurefunc(QuantumState state1, QuantumState state2, int i, int j, int k, int l);         // basic function for next two functions

    // s-orbital functions

    //chem. potentials
    // DEPRECATED
    void add_mu(int st, RealType mu);                        //adds chem. potential on multiorbital
    // DEPRECATED
    void add_mus(int st, RealType mus);                        //adds chem. potential on s-orbitals

    //functions describe hoppings
    // DEPRECATED
    int checkhop(long int state1, long int state2, int i, int j);        //check probability hopping between state1 and state2
    // DEPRECATED
    void add_hopping(RealMatrixType& HoppingMatrix);            //function adds to Hamilt hopping electron from "i"
    // DEPRECATED
    void add_hopping(int i,int j, RealType t);            //function adds to Hamilt hopping electron from "i"

public:

    HamiltonianPart(IndexClassification &F, StatesClassification &S, const QuantumNumbers& QN);

    void prepare(void);
    void diagonalize(void);
    bool reduce(RealType ActualCutoff);

    InnerQuantumState getSize(void) const;

    RealType getMatrixElement(InnerQuantumState m, InnerQuantumState n) const; //return H(m,n)
    RealType getEigenValue(InnerQuantumState state) const; //return V(m)
    RealType getMinimumEigenvalue() const;        //!<Return the lowest Eigenvalue of the current part;
    RealVectorType getEigenState(InnerQuantumState state) const;

    QuantumNumbers getQuantumNumbers() const;                //return quantum numbers of current Hamiltonian part
    BlockNumber getBlockNumber() const;

    // DEPRECATED
    void print_to_screen();            //print to screen part of hamiltonian
    
    void save(H5::CommonFG* RootGroup) const;
    void load(const H5::CommonFG* RootGroup);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIANPART_H
