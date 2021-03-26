/** \file include/pomerol/FieldOperatorPart.h
** \brief Declaration of FieldOperatorPart, CreationOperatorPart and AnnihilationOperatorPart classes.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATORPART_H
#define __INCLUDE_FIELDOPERATORPART_H

#include "Misc.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "OperatorPresets.h"

namespace Pomerol{

template<bool Complex> class FieldOperator;
template<bool Complex> class FieldOperatorContainer;

/** This class is an abstract implementation of the electronic creation/annihilation operators, which acts in the eigenbasis of the Hamiltonian
 * between it's certain blocks.
 * Rotation to the basis is done in the following way:
 * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
 * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
 * The actual creation and annihilation operators are inherited.
 */
template<bool Complex = false>
class FieldOperatorPart : public ComputableObject {
    friend class FieldOperator<Complex>;
    friend class FieldOperatorContainer<Complex>;
public:

    /** A reference to the IndexClassification object. */
    const IndexClassification<Complex> &IndexInfo;
    /** A reference to the StateClassification object. */
    const StatesClassification<Complex> &S;
    /** A reference to the HamiltonianPart on the right hand side. */
    const HamiltonianPart<Complex> &HFrom;
    /** A reference to the HamiltonianPart on the left hand side. */
    const HamiltonianPart<Complex> &HTo;
protected:
    /** A pointer to the Operator object ( Pomerol::OperatorPresets::C or Cdag ). */
    Operator<Complex> *O;

    /** Index of the field operator. */
    ParticleIndex PIndex;
    /** Storage of the matrix elements of the operator. Row ordered sparse matrix. */
    RowMajorMatrixType<Complex> elementsRowMajor;
    /** Copy of the Storage of the matrix elements of the operator. Column ordered sparse matrix. */
    ColMajorMatrixType<Complex> elementsColMajor;
    /** The tolerance with which the matrix elements are evaluated. */
    //static
    const RealType MatrixElementTolerance; //1e-8 by default
    /** Make this class purely abstract. */
    virtual void do_nothing(void) = 0;

public:

    /** Constructor.
     * \param[in] IndexInfo A const reference to the IndexClassification object.
     * \param[in] S A const reference to the StateClassification object.
     * \param[in] HFrom A const reference to the HamiltonianPart on the right hand side.
     * \param[in] HTo A const reference to the HamiltonianPart on the left hand side.
     * \param[in] PIndex Index of the field operator.
     */
    FieldOperatorPart(const IndexClassification<Complex> &IndexInfo,
                      const StatesClassification<Complex> &S,
                      const HamiltonianPart<Complex> &HFrom,
                      const HamiltonianPart<Complex> &HTo,
                      ParticleIndex PIndex);

    /** Compute all the matrix elements. Changes the Status of the object to Computed. */
    void compute();
    /** Print all matrix elements of the operator to screen. */
    void print_to_screen() const;

    /** Returns the row ordered sparse matrix of matrix elements. */
    const RowMajorMatrixType<Complex>& getRowMajorValue(void) const;
    /** Returns the column ordered sparse matrix of matrix elements. */
    const ColMajorMatrixType<Complex>& getColMajorValue(void) const;
    /** Returns the right hand side index. */
    BlockNumber getRightIndex(void) const;
    /** Returns the left hand side index. */
    BlockNumber getLeftIndex(void) const;

    #ifdef ENABLE_SAVE_PLAINTEXT
    /** Save the data to the file.
     * \param[in] path Path to the file.
     */
    bool savetxt(const boost::filesystem::path &path);
    #endif
};

// Forward declarations
template<bool Complex> class AnnihilationOperatorPart;
template<bool Complex> class CreationOperatorPart;
template<bool Complex> class QuadraticOperatorPart;

template<bool Complex> class AnnihilationOperator;
template<bool Complex> class CreationOperator;
template<bool Complex> class QuadraticOperator;

/** This class is inherited from FieldOperatorPart and is a part of electronic annihilation operator in the eigenbasis of the Hamiltonian between it's two blocks. */
template<bool Complex = false>
class AnnihilationOperatorPart : public FieldOperatorPart<Complex>
{
    friend class CreationOperatorPart<Complex>;
    friend class CreationOperator<Complex>;
    friend class AnnihilationOperator<Complex>;
    /** Does nothing. Private. */
    void do_nothing(){};
public :
    /** Constructor. Look FieldOperatorPart::FieldOperatorPart. */
    AnnihilationOperatorPart(const IndexClassification<Complex> &IndexInfo,
                             const StatesClassification<Complex> &S,
                             const HamiltonianPart<Complex> &HFrom,
                             const HamiltonianPart<Complex> &HTo,
                             ParticleIndex PIndex);
    /** Construct the CreationOperatorPart from the class (transpose it). */
    const CreationOperatorPart<Complex>& transpose(void) const;
};

/** This class is inherited from FieldOperatorPart and is a part of electronic creation operator in the eigenbasis of the Hamiltonian between it's two blocks. */
template<bool Complex = false>
class CreationOperatorPart : public FieldOperatorPart<Complex>
{
    friend class AnnihilationOperatorPart<Complex>;
    friend class AnnihilationOperator<Complex>;
    friend class CreationOperator<Complex>;
    /** Does nothing. Private. */
    void do_nothing(){};
public :
    /** Constructor. Look FieldOperatorPart::FieldOperatorPart. */
    CreationOperatorPart(const IndexClassification<Complex> &IndexInfo,
                         const StatesClassification<Complex> &S,
                         const HamiltonianPart<Complex> &HFrom,
                         const HamiltonianPart<Complex> &HTo,
                         ParticleIndex PIndex);
    /** Construct the AnnihilationOperatorPart from the class (transpose it). */
    const AnnihilationOperatorPart<Complex>& transpose(void) const;
};

template<bool Complex = false>
class QuadraticOperatorPart : public FieldOperatorPart<Complex>
{
    friend class QuadraticOperator<Complex>;
    /** Does nothing. Private. */
    void do_nothing(){};
protected:
    /** Indices of the operator. Used instead of FieldOperator::Index */
    ParticleIndex Index1, Index2;
public:
    QuadraticOperatorPart(const IndexClassification<Complex> &IndexInfo,
                          const StatesClassification<Complex> &S,
                          const HamiltonianPart<Complex> &HFrom,
                          const HamiltonianPart<Complex> &HTo,
                          ParticleIndex PIndex1, ParticleIndex PIndex2);
};

// External templates: Real case

extern template class FieldOperatorPart<false>;
extern template class AnnihilationOperatorPart<false>;
extern template class CreationOperatorPart<false>;
extern template class QuadraticOperatorPart<false>;

// External templates: Complex case

extern template class FieldOperatorPart<true>;
extern template class AnnihilationOperatorPart<true>;
extern template class CreationOperatorPart<true>;
extern template class QuadraticOperatorPart<true>;

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H
