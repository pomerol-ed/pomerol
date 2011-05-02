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


#ifndef __INCLUDE_STATESCLASSIFICATION_H
#define __INCLUDE_STATESCLASSIFICATION_H

#include"Misc.h"
#include"HDF5Storage.h"
#include"IndexClassification.h"

/** This class represents a number of a current block of QuantumStates which have
 * the same quantum numbers. If such block cannot exist it's value is assigned to an ERROR_BLOCK_NUMBER
 * The classification of blocks is now done by StatesClassification class
 */
struct BlockNumber {
    
    /** The number of current block */
    int number;
    
    /** Empty constructor */
    BlockNumber(){};

    /** Copy constructor from int number
     * \param[in] number A number to copy from 
     */
    BlockNumber(int number):number(number){};

    /** Type conversion to integer */
    operator int() const {return number;}

    /** post-incremet operator */
    BlockNumber& operator ++(int unused){number++; return *this;}

    /** Returns true if such block exists */
    bool isCorrect(){return number >= 0;}
    /** Operator < */
    bool operator<(const BlockNumber& rhs) const ;
    /** Operator > */
    bool operator==(const BlockNumber& rhs) const;
};

/** All blocks with this number are treated as nonexistent */
const BlockNumber ERROR_BLOCK_NUMBER = -1;

/** This struct is a set of QuantumNumbers */
// TODO: remove explicit dependence from Lz, N_up, N_down
struct QuantumNumbers : public std::vector<short> {

    //short Lz; 	[0]
    //short N_up;	[1]
    //short N_down;	[2]

    QuantumNumbers(int LZ, int N_UP, int N_DOWN);
    /** Empty constructor */
    QuantumNumbers();
    /** streaming << operator */
    friend std::ostream& operator<<(std::ostream& output, const QuantumNumbers& out);
    /** == operator 
     * \param[in] rhs A right hand side of an operator 
     */
    bool operator==(const QuantumNumbers &rhs)const ;
    /* < operator 
     * \param[in] rhs A right hand side of an operator 
     */
    bool operator<(const QuantumNumbers &rhs) const ;

    /** Save to an output hdf5 file
     * \param[in] RootGroup Group ih h5 file under which to save the file
     */
    //void save(H5::CommonFG* RootGroup) const;
    /** Load a value from hdf5 file
     * \param[in] RootGroup Group ih h5 file under which QuantumNumbers is stored
     */
    //void load(const H5::CommonFG* RootGroup);
};

/** All blocks with current QuantumNumber are treated as non-existent */
const QuantumNumbers ERROR_QUANTUM_NUMBERS = QuantumNumbers(0,-1,-1);

/** A typedef for a quantum state in the Hilbert space */
typedef unsigned long int QuantumState;
/** A typedef for a number of a QuantumState in a current block. In case 
 * of no symmetry and only 1 block InnerQuantum state would be equivalent to QuantumState */
typedef unsigned long int InnerQuantumState;


/** This class handles all information about Quantum States in a Fock basis. 
 * It defines all blocks of QuantumStates and classifies the Quantum States
 * It is called prior to Hamiltonian, since it relies only IndexInfo
 */
class StatesClassification {
    
    /** Total number of states = 2^(IndexInfo.size()) */
    QuantumState N_state;            

    /** Total number of modes of the system. Equal to IndexInfo.size() */
    int N_bit;                //number bit of states
    // TODO : remove N_bit_m
    int N_bit_m;                //2*(2*L(orbital) +1)
    // TODO : remove ***st
    std::vector<QuantumState> *** st;        
    // TODO : remove size
    int size;                //number of classificated vectors
    /* A map between all BlockNumbers and QuantumNumbers */
    std::map<BlockNumber,QuantumNumbers> BlockToQuantum;
    /* A map between all QuantumNumbers and BlockNumbers */
    std::map<QuantumNumbers,BlockNumber> QuantumToBlock;

    // TODO : remove maximumBlockNumber_
    BlockNumber maximumBlockNumber_; 
    // TODO : remove BLOCKNUMBERLIMIT
    unsigned int BLOCKNUMBERLIMIT;

    /** A reference to an IndexClassification object */
    IndexClassification &IndexInfo;

    // TODO: remove putstates()
    void putstates();            //function gets all clasificated states with Lz = "Lz", N_up = "N_up", N_down = "N_down"
public:        
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     */
    StatesClassification(IndexClassification& IndexInfo):IndexInfo(IndexInfo) {};

    /** Perform a classification of all QuantumStates */
    void compute();                     

    const int N_b_m();                            //return N_bit_m
    const int L();                                //return value of orbital moment
    const QuantumState N_st();                            //return N_state

    const std::vector<QuantumState>& clstates( QuantumNumbers in );                //return st[in.Lz][in.N_up][in.N_down]
    const QuantumState cst( QuantumNumbers in, int m);                    //return st[in.Lz][in.N_up][in.N_down][m]
    const InnerQuantumState getInnerState( QuantumState state);                    //finds number of state @state in corresponding block 

    BlockNumber getBlockNumber(QuantumNumbers in);            //returns a number of Block which corresponds to given Quantum Numbers
    QuantumNumbers getBlockInfo(BlockNumber in);        //return Lz,N_up,N_down of number num
    BlockNumber NumberOfBlocks();                        //return amount of non-trivial hamiltonian blocks
    QuantumNumbers getStateInfo(QuantumState in);        //return Lz,N_up,N_down of number num
    BlockNumber getBlockNumber(QuantumState in);            //returns a number of Block which corresponds to given Quantum Numbers

    inline int n_i(long int state, int i) { return ((state&(1<<i))>>i) ;};
    void getSiteInfo(int bit, int& lz, int& spin);

    bool checkQuantumNumbers(QuantumNumbers in);
};

#endif // endif :: #ifndef __INCLUDE_STATESCLASSIFICATION_H
