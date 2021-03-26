/** \file include/pomerol/FieldOperator.h
** \brief Declaration of field operators : creation and annihilation operators.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATOR_H
#define __INCLUDE_FIELDOPERATOR_H

#include<boost/bimap.hpp>

#include <mpi.h>

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
template<bool Complex = false>
class FieldOperator : public ComputableObject
{
public:

    typedef boost::bimaps::bimap<
        boost::bimaps::set_of<BlockNumber>,
        boost::bimaps::set_of<BlockNumber>
    > BlocksBimap;
    typedef BlocksBimap::value_type BlockMapping;

    using PartT = FieldOperatorPart<Complex>;

protected:
    /** A reference to a IndexClassification object */
    const IndexClassification<Complex> &IndexInfo;
    /** A reference to a StatesClassification object */
    const StatesClassification<Complex> &S;
    /** A reference to a Hamiltonian object */
    const Hamiltonian<Complex> &H;
    /** A reference an Operator object (OperatorPresets::C or Cdag). */
    const Operator<Complex> *O;

    /** An index of the operator */
    ParticleIndex Index;
    /** A vector of parts */
    std::vector<PartT*> parts;
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
    virtual QuantumNumbers<Complex> mapsTo(const QuantumNumbers<Complex>& in) const;

public:
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    FieldOperator(const IndexClassification<Complex> &IndexInfo,
                  const StatesClassification<Complex> &S,
                  const Hamiltonian<Complex> &H,
                  ParticleIndex Index);

    /** Returns a FieldOperatorPart based on its left BlockNumber */
    PartT& getPartFromLeftIndex(BlockNumber in) const;
    /** Returns a FieldOperatorPart based on its left QuantumNumbers */
    PartT& getPartFromLeftIndex(const QuantumNumbers<Complex>& in) const;
    /** Returns a FieldOperatorPart based on its right BlockNumber */
    PartT& getPartFromRightIndex(BlockNumber out) const;
    /** Returns a FieldOperatorPart based on its right QuantumNumbers */
    PartT& getPartFromRightIndex(const QuantumNumbers<Complex>& out) const;
    /** Returns a left BlockNumber for a given right BlockNumber */
    BlockNumber getLeftIndex(BlockNumber RightIndex) const;
    /** Returns a right BlockNumber for a given left BlockNumber */
    BlockNumber getRightIndex(BlockNumber LeftIndex) const;
    /** Returns a reference to BlockMapping */
    BlocksBimap const& getBlockMapping() const;

    /** Returns a vector of all underlying parts. */
    const std::vector<PartT*>& getParts();

    /** Returns acting ParticleIndex of current operator */
    ParticleIndex getIndex(void) const;
    /** Virtual method for assigning world-lines */
    virtual void prepare(void) = 0;
    /** Computes all world-lines */
    void compute(const MPI_Comm& comm = MPI_COMM_WORLD);
};

/** A creation operator in the eigenspace of a Hamiltonian */
template<bool Complex> class CreationOperator;
/** An annihilation operator in the eigenspace of a Hamiltonian */
template<bool Complex> class AnnihilationOperator;
/** A quadratic operator, c_1^+ c_2, in the eigenspace of a Hamiltonian */
template<bool Complex> class QuadraticOperator;

template<bool Complex = false>
class CreationOperator : public FieldOperator<Complex>
{
    friend class AnnihilationOperator<Complex>;
    friend class FieldOperatorContainer<Complex>;
    friend class QuadraticOperator<Complex>;
public:

    using Base = FieldOperator<Complex>;

    /* Returns hermitian conjugate of current operator */
    AnnihilationOperator<Complex>& transpose(void);
    void prepare();

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    CreationOperator(const IndexClassification<Complex> &IndexInfo,
                     const StatesClassification<Complex> &S,
                     const Hamiltonian<Complex> &H,
                     ParticleIndex Index);
};

template<bool Complex = false>
class AnnihilationOperator : public FieldOperator<Complex>
{
    friend class CreationOperator<Complex>;
    friend class FieldOperatorContainer<Complex>;
    friend class QuadraticOperator<Complex>;
public:

    using Base = FieldOperator<Complex>;

    /* Returns hermitian conjugate of current operator */
    CreationOperator<Complex>& transpose(void);

    void prepare();

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index An index of an operator
     */
    AnnihilationOperator(const IndexClassification<Complex> &IndexInfo,
                         const StatesClassification<Complex> &S,
                         const Hamiltonian<Complex> &H,
                         ParticleIndex Index);
};

template<bool Complex = false>
class QuadraticOperator : public FieldOperator<Complex>
{
protected:
    /** Indices of the operator. Used instead of FieldOperator::Index */
    ParticleIndex Index1, Index2;

public:

    using Base = FieldOperator<Complex>;

    void prepare();

    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     * \param[in] S A reference to a StatesClassification object
     * \param[in] H A reference to a Hamiltonian object
     * \param[in] Index1 An index of a creation operator
     * \param[in] Index2 An index of an annihilation operator
     */
    QuadraticOperator(const IndexClassification<Complex> &IndexInfo,
                      const StatesClassification<Complex> &S,
                      const Hamiltonian<Complex> &H,
                      ParticleIndex Index1, ParticleIndex Index2);
};

// External templates: Real case

extern template class FieldOperator<false>;
extern template class AnnihilationOperator<false>;
extern template class CreationOperator<false>;
extern template class QuadraticOperator<false>;

// External templates: Complex case

extern template class FieldOperator<true>;
extern template class AnnihilationOperator<true>;
extern template class CreationOperator<true>;
extern template class QuadraticOperator<true>;

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATOR_H
