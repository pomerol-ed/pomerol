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


#ifndef __INCLUDE_HAMILTONIANPART_H
#define __INCLUDE_HAMILTONIANPART_H
#include "Misc.h"
#include "HDF5Storage.h"
#include "ComputableObject.h"
#include "IndexClassification.h"
#include "StatesClassification.h"

namespace Pomerol{

class HamiltonianPart : public ComputableObject, public HDF5Storable {

    IndexClassification &IndexInfo;
    StatesClassification &S;

    QuantumNumbers hpart_id;

    RealMatrixType H;                //part of Hamiltonian
    RealVectorType V;                //vector of Eigen Values

    void add_nTerm(InnerQuantumState st, nTerm *N);

    void add_nnTerm(InnerQuantumState st, nnTerm *T);
    void add_spinflipTerm(InnerQuantumState st, spinflipTerm *T);

    int measurefunc(QuantumState state1, QuantumState state2, int i, int j, int k, int l);         // basic function for next two functions

    // s-orbital functions

    //chem. potentials

    void add_mu(int st, RealType mu);                        //adds chem. potential on multiorbital
    void add_mus(int st, RealType mus);                        //adds chem. potential on s-orbitals

    //functions describe hoppings

    int checkhop(long int state1, long int state2, int i, int j);        //check probability hopping between state1 and state2
    void add_hopping(RealMatrixType& HoppingMatrix);            //function adds to Hamilt hopping electron from "i"
    void add_hopping(int i,int j, RealType t);            //function adds to Hamilt hopping electron from "i"

public:

    HamiltonianPart(IndexClassification &F, StatesClassification &S, QuantumNumbers id);

    void enter();

    InnerQuantumState size(void);
    RealType reH(int m, int n);        //return H(m,n)
    RealType reV(int m);            //return V(m)

    void diagonalization();                //method of process diagonalization
    QuantumNumbers id();                //return id of current hpart
    BlockNumber getId();
    RealType getMinimumEigenvalue();        //!<Return the lowest Eigenvalue of the current part;
    
    bool reduce(RealType ActualCutoff);
    void print_to_screen();            //print to screen part of hamiltonian
    RealVectorType getEigenState(InnerQuantumState m);

    void save(H5::CommonFG* RootGroup) const;
    void load(const H5::CommonFG* RootGroup);

};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIANPART_H
