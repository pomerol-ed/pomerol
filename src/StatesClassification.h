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


/** \file src/StatesClassification.h
** \brief Declaration of BlockNumber, QuantumNumbers and StatesClassification classes
**
** \author Andrey Antipov (antipov@ct-qmc.org)
** \author Mikhail Aleynikov (alejnikov89@mail.ru)
** \author Igor Krivenko (igor@shg.ru)
*/


#ifndef __INCLUDE_STATESCLASSIFICATION_H
#define __INCLUDE_STATESCLASSIFICATION_H

#include"Misc.h"
#include"HDF5Storage.h"
#include"IndexClassification.h"

namespace Pomerol{

/** This class represents a number of a current block of QuantumStates which have
 * the same quantum numbers. If such block can't exist ( can't be created by anyone else )
 * it's value is assigned to an ERROR_BLOCK_NUMBER.
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

    /** post-increment operator */
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

/** This struct is a set of quantum numbers, available in current model */
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

    /** get total number of Quantum States ( 2^IndexInfo.size() ) */
    const QuantumState getNumberOfStates();                           

    /** get a vector of all QuantumStates with a given set of QuantumNumbers
     * \param[in] in A set of quantum numbers to get a vector of QuantumStates 
     */
    const std::vector<QuantumState>& getQuantumStates( QuantumNumbers in );   

    /** get a QuantumState, corresponding to an internal InnerQuantumState
     * \param[in] QuantumNumbers of block in which the InnerQuantumState is located
     * \param[in] m InnerQuantumState for which the correspondence is required
     */
    const QuantumState getQuantumState( QuantumNumbers in, int m);
    /** get InnerQuantumState of a given QuantumState. Since QuantumState is associated with
     * the Block number no explicit BlockNumber or QuantumNumbers is required 
     * \param[in] state QuantumState for which the correspondence is required
     */
    const InnerQuantumState getInnerState( QuantumState state);   

    /** Returns a number of Block which corresponds to given Quantum Numbers 
     * \param[in] in A set of QuantumNumbers to find corresponding BlockNumber
     */
    BlockNumber getBlockNumber(QuantumNumbers in);            

    /** Returns QuantumNumbers for a given BlockNumber
     * \param[in] in A BlockNumber to find a set of corresponding QuantumNumbers
     */
    QuantumNumbers getBlockInfo(BlockNumber in);        
    /** Returns total amount of non-vanishing blocks */
    BlockNumber NumberOfBlocks();                        

    /** Returns QuantumNumbers of a given QuantumState 
     * \param[in] in A QuantumState for which the QuantumNumbers are requested
     */
    QuantumNumbers getStateInfo(QuantumState in);        
    /** Returns BlockNumber of a given QuantumState 
     * \param[in] in A QuantumState for which the BlockNumber is requested
     */
    BlockNumber getBlockNumber(QuantumState in);       

    /** Returns the value of <QuantumState|\hat n_i|QuantumState> */
    inline int n_i(QuantumState state, ParticleIndex i) { return ((state&(1<<i))>>i) ;};
    // TODO: rewrite this method
    void getSiteInfo(int bit, int& lz, int& spin);

    /** Checks that a block with a given QuantumNumbers does not vanish 
     * \param[in] in A set of QuantumNumbers to check
     */
    bool checkQuantumNumbers(QuantumNumbers in);
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_STATESCLASSIFICATION_H
