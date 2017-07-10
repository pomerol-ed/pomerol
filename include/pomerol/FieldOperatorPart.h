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

/** This class is an abstract implementation of the electronic creation/annihilation operators, which acts in the eigenbasis of the Hamiltonian
 * between it's certain blocks.
 * Rotation to the basis is done in the following way:
 * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
 * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
 * The actual creation and annihilation operators are inherited.
 */
class FieldOperatorPart : public ComputableObject {
    friend class FieldOperatorContainer;
    friend class FieldOperator;
public:
    /** A reference to the IndexClassification object. */
    const IndexClassification &IndexInfo;
    /** A reference to the StateClassification object. */
    const StatesClassification &S;
    /** A reference to the HamiltonianPart on the right hand side. */
    const HamiltonianPart &HFrom;
    /** A reference to the HamiltonianPart on the left hand side. */
    const HamiltonianPart &HTo;
protected:
    /** A pointer to the Operator object ( Pomerol::OperatorPresets::C or Cdag ). */
    Operator *O;

    /** Index of the field operator. */
    ParticleIndex PIndex;
    /** Storage of the matrix elements of the operator. Row ordered sparse matrix. */
    RowMajorMatrixType elementsRowMajor;
    /** Copy of the Storage of the matrix elements of the operator. Column ordered sparse matrix. */
    ColMajorMatrixType elementsColMajor;
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
    FieldOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);

    /** Compute all the matrix elements. Changes the Status of the object to Computed. */
    void compute();
    /** Print all matrix elements of the operator to screen. */
    void print_to_screen() const;

    /** Returns the row ordered sparse matrix of matrix elements. */
    const RowMajorMatrixType& getRowMajorValue(void) const;
    /** Returns the column ordered sparse matrix of matrix elements. */
    const ColMajorMatrixType& getColMajorValue(void) const;
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
class AnnihilationOperatorPart;
class CreationOperatorPart;

/** This class is inheried from FieldOperatorPart and is a part of electronic annihilation operator in the eigenbasis of the Hamiltonian between it's two blocks. */
class AnnihilationOperatorPart : public FieldOperatorPart
{
    friend class CreationOperatorPart;
    friend class CreationOperator;
    friend class AnnihilationOperator;
    /** Does nothing. Private. */
    void do_nothing(){};
public :
    /** Constructor. Look FieldOperatorPart::FieldOperatorPart. */
    AnnihilationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);
    /** Construct the CreationOperatorPart from the class (transpose it). */
    const CreationOperatorPart& transpose(void) const;
};

/** This class is inheried from FieldOperatorPart and is a part of electronic creation operator in the eigenbasis of the Hamiltonian between it's two blocks. */
class CreationOperatorPart : public FieldOperatorPart
{
    friend class AnnihilationOperatorPart;
    friend class AnnihilationOperator;
    friend class CreationOperator;
    /** Does nothing. Private. */
    void do_nothing(){};
public :
    /** Constructor. Look FieldOperatorPart::FieldOperatorPart. */
    CreationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex);
    /** Construct the AnnihilationOperatorPart from the class (transpose it). */
    const AnnihilationOperatorPart& transpose(void) const;
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H
