//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file src/FieldOperator.h
** \brief Declaration of field operators : creation and annihilation operators.
** 
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATOR_H
#define __INCLUDE_FIELDOPERATOR_H

#include<boost/bimap.hpp>

#include"Misc.h"
#include"StatesClassification.h"
#include"Hamiltonian.h" 
#include"FieldOperatorPart.h"

namespace Pomerol{

/** \typedef 
 * A pair of left and right indices of a part in a Field Operator. Each part is a non-vanishing worldline in an operator
 */
typedef std::pair<BlockNumber,BlockNumber> BlockMapping;

/** This class is a parent class for creation/annihilation operators which act
 * on all blocks of quantum states */ 
class FieldOperator : public ComputableObject 
{
public:
    
    typedef boost::bimaps::bimap<
        boost::bimaps::set_of<BlockNumber>,
        boost::bimaps::set_of<BlockNumber>
    > BlocksBimap;
    typedef BlocksBimap::value_type BlockMapping;
    
protected:
    /** A reference to a IndexClassification object */
    const IndexClassification &IndexInfo;
    /** A reference to a StatesClassification object */
    const StatesClassification &S;
    /** A reference to a Hamiltonian object */
    const Hamiltonian &H;
    /** A reference an Operator object (OperatorPresets::C or Cdag). */
    const Operator *O;

    /** An index of the operator */
    ParticleIndex Index;
    /** A vector of parts */
    std::vector<FieldOperatorPart*> parts;
    /** A map between non-vanishing parts (internal numbering) and their R.H.S. BlockNumbers  */
    std::map<size_t,BlockNumber> mapPartsFromRight;
    /** A map between non-vanishing parts (internal numbering) and their L.H.S. BlockNumbers  */
    std::map<size_t,BlockNumber> mapPartsFromLeft;
        
    BlocksBimap LeftRightBlocks;

    /** Return the resulting BlockNumber of states obtained by this operator, acting on states from another block. 
     * If no BlockNumber found returns ERROR_BLOCK_NUMBER.
     * \param[in] RightIndex The BlockNumber of states on right hand side of the FieldOperator.
     */
    virtual BlockNumber mapsTo(BlockNumber RightIndex) const;

    /** Return the resulting QuantumNumbers of states obtained by this operator, acting on states from another block. 
     * If no QuantumNumbers found throws an exception.
     * \param[in] RightIndex The BlockNumber of states on right hand side of the FieldOperator.
     */
    virtual QuantumNumbers mapsTo(const QuantumNumbers& in) const;

public:
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    FieldOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index);

    /** Returns a FieldOperatorPart based on its left BlockNumber */
    const FieldOperatorPart& getPartFromLeftIndex(BlockNumber in) const;
    /** Returns a FieldOperatorPart based on its left QuantumNumbers */
    const FieldOperatorPart& getPartFromLeftIndex(const QuantumNumbers& in) const;
    /** Returns a FieldOperatorPart based on its right BlockNumber */
    const FieldOperatorPart& getPartFromRightIndex(BlockNumber out) const;
    /** Returns a FieldOperatorPart based on its right QuantumNumbers */
    const FieldOperatorPart& getPartFromRightIndex(const QuantumNumbers& out) const;
    /** Returns a left BlockNumber for a given right BlockNumber */
    BlockNumber getLeftIndex(BlockNumber RightIndex) const;
    /** Returns a right BlockNumber for a given left BlockNumber */
    BlockNumber getRightIndex(BlockNumber LeftIndex) const;
    /** Returns a reference to BlockMapping */
    BlocksBimap const& getBlockMapping() const;

    /** Returns acting ParticleIndex of current operator */
    ParticleIndex getIndex(void) const;
    /** Virtual method for assigning world-lines */
    virtual void prepare(void) = 0;
    /** Computes all world-lines */
    void compute(void);
};

/** A creation operator in the eigenspace of a Hamiltonian */
class CreationOperator;
/** An annihilation operator in the eigenspace of a Hamiltonian */
class AnnihilationOperator;

class CreationOperator : public FieldOperator
{
    friend class AnnihilationOperator;
public:
    /* Returns hermitian conjugate of current operator */
    AnnihilationOperator& transpose(void);
    void prepare();

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    CreationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index);
};

class AnnihilationOperator : public FieldOperator
{
    friend class CreationOperator;
public:
    /* Returns hermitian conjugate of current operator */
    CreationOperator& transpose(void);

    void prepare();

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    AnnihilationOperator(const IndexClassification &IndexInfo, const StatesClassification &S, const Hamiltonian &H, ParticleIndex Index);
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATOR_H
