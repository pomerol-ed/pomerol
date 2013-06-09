//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file Symmetrizer.h
**  \brief Declaration of the Symmetrizer class - a class to store and get the information about the symmetries of the system.
** 
**  \author    Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_SYMMETRIZER_H
#define __INCLUDE_SYMMETRIZER_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "IndexHamiltonian.h"
#include "ComputableObject.h"
#include <boost/functional/hash.hpp>
#include <set>

namespace Pomerol{


/** This class stores the information about operations, which commute with the Hamiltonian.
 * It tries to find Lattice symmetries, checks for some common symmetries 
 * and also checks given symmetries. */
class Symmetrizer : public ComputableObject
{
//typedef void (Symmetrizer::*OperatorPtr)(Operator *, Json::Value&);
public:
    struct IndexPermutation;
    /** This generates a trivial combination of 0123...N-1 indices. */
    static const DynamicIndexCombination& generateTrivialCombination(ParticleIndex N);
    /** This class represents a set of conserved quantum numbers. It is accompanied by a list of operators which should act on the state. */
    struct QuantumNumbers;
    /** Statuses of the object */
    enum {Constructed,Computed};
private:
    /** A link to an IndexClassification object. */ 
    const IndexClassification &IndexInfo;
    /** A link to an IndexHamiltonian object. */
    const IndexHamiltonian &Storage;

    /** Total amount of indices in the system. */
    const ParticleIndex IndexSize;

    /** A list of equivalent lattice sites permutations. */
    std::list<IndexPermutation*> Permutations;
    /** Total amount of symmetries found. */
    int NSymmetries;
    /** A vector of operators that commute with the Hamiltonian. */
    std::vector<boost::shared_ptr<Operator> > Operations;


    /** This method finds all possible symmetry operations. */ /** lattice permutation operators, that commute with the hamiltonian. */
    //void findLatticeSymmetry();
public:
    Symmetrizer(IndexClassification &IndexInfo, IndexHamiltonian &Storage);
    /** This method checks several possible symmetry operations to split the Hamiltonian into blocks. */
    void compute(bool ignore_symmetries = false);

    bool checkSymmetry(const Operator &in);
    /** Get a vector of operators that commute with the Hamiltonian. */
    const std::vector<boost::shared_ptr<Operator> >& getOperations() const;
    /** Get a sample QuantumNumbers. Their amount is set. */
    QuantumNumbers getQuantumNumbers() const;
};

/** A combination of indices to which a permutation commutes with a Hamiltonian. 
 * Only the combinations with all different indices can be accepted. Also, since it
 * represents a permutation of indices a check for the irreducibility of the permutation is done, 
 * i.e. it is checked that a permutation can not be split in at least two others. 
 * It is assumed that a trivial identity permutation makes an error. */
struct Symmetrizer::IndexPermutation 
{
private:
    std::vector<DynamicIndexCombination*> Combinations;
    /** This defines how many permutations should be done to form a closed loop. */
    unsigned int CycleLength;
    /** Calculates CycleLength. */
    void calculateCycleLength();
    /** Checks that all elements of permutation are different and belong to 0..N-1 interval. */
    bool checkConsistency(const DynamicIndexCombination &in);
    /** Checks that the permutation can not be splitted. */
    bool checkIrreducibility(const DynamicIndexCombination &in); 
    /** Total amount of indices. */
    const ParticleIndex N;
public:
    /** Constructor. 
     * \param[in] @in Equivalent index combination. 
     */
    IndexPermutation(const DynamicIndexCombination &in);

    /** Returns the permutation for a given number in cycle
     * \param[in] cycle Defines which permutation of indices to return. By default returns the first one.
     */
    const DynamicIndexCombination& getIndices( unsigned int cycle = 1) const;

    /** Returns the length of the permutation cycle. */
    const unsigned int getCycleLength() const;
    /** Exception - equal indices. */
    class exEqualIndices : public std::exception { virtual const char* what() const throw(); };

};

/** This class represents a set of quantum numbers obtained by the Symmetrizer. */
struct Symmetrizer::QuantumNumbers { 
friend class Symmetrizer;
private:
    /** Total number of quantum numbers. */
    int amount;
    /** A vector of numbers. For now set as MelemType. */
    std::vector<MelemType> numbers;
    /** Private constuctor - can be called only inside Symmetrizer. */
    QuantumNumbers(int amount);
    /** A hash of the quantum numbers - used for comparison. */
    std::size_t NumbersHash;
    /** Hash generator. */
    boost::hash<std::vector<MelemType> > numbers_hash_generator;
public:
    /** Set a quantum number at the given position to a value
     * \param[in] pos Position of QuantumNumber, e.g. the number of operation in Symmetrizer::Operations.
     * \param[in] val Value of QuantumNumber.
     */
    bool set ( int pos, MelemType val );

    /* Comparison operators. */
    bool operator< (const QuantumNumbers& rhs) const ;
    bool operator== (const QuantumNumbers& rhs) const ;
    bool operator!= (const QuantumNumbers& rhs) const ;
    /** Output to external stream */
    friend std::ostream& operator<<(std::ostream& output, const QuantumNumbers& out);
    /** Exception for bad quantum numbers. */
    class exWrongNumbers : public std::exception { virtual const char* what() const throw() { return "Wrong QuantumNumbers."; } };
};

/** A typedef for QuantumNumbers since it is often used. */ 
typedef Symmetrizer::QuantumNumbers QuantumNumbers;

}; // end of namespace Pomerol

#endif //  endif :: #ifndef __INCLUDE_SYMMETRIZER_H
