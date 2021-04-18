/** \file include/pomerol/FieldOperatorPart.h
** \brief Declaration of FieldOperatorPart, CreationOperatorPart and AnnihilationOperatorPart classes.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_FIELDOPERATORPART_H
#define __INCLUDE_FIELDOPERATORPART_H

#include "Misc.h"
#include "HilbertSpace.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
//#include "OperatorPresets.h"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <memory>

namespace Pomerol {

/** This class is an abstract implementation of the electronic creation/annihilation operators, which acts in the eigenbasis of the Hamiltonian
 * between it's certain blocks.
 * Rotation to the basis is done in the following way:
 * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
 * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
 * The actual creation and annihilation operators are inherited.
 */
class MonomialOperatorPart : public ComputableObject {
    friend class FieldOperatorContainer;
    friend class FieldOperator;

protected:
    template<bool C>
    using LOperatorType = libcommute::loperator<MelemType<C>, libcommute::fermion>;

private:
    bool Complex;

    const void* MOp;

    /** A reference to the StateClassification object. */
    const StatesClassification &S;
    /** A reference to the HamiltonianPart on the right hand side. */
    const HamiltonianPart &HFrom;
    /** A reference to the HamiltonianPart on the left hand side. */
    const HamiltonianPart &HTo;
protected:

    /** Storage of the matrix elements of the operator. Row ordered sparse matrix. */
    std::shared_ptr<void> elementsRowMajor;
    /** Copy of the Storage of the matrix elements of the operator. Column ordered sparse matrix. */
    std::shared_ptr<void> elementsColMajor;
    /** The tolerance with which the matrix elements are evaluated. */
    //static
    const RealType MatrixElementTolerance = 1e-8;

public:

    /** Constructor.
     * \param[in] IndexInfo A const reference to the IndexClassification object.
     * \param[in] S A const reference to the StateClassification object.
     * \param[in] HFrom A const reference to the HamiltonianPart on the right hand side.
     * \param[in] HTo A const reference to the HamiltonianPart on the left hand side.
     * \param[in] PIndex Index of the field operator.
     */
    template<typename ScalarType>
    MonomialOperatorPart(const libcommute::loperator<ScalarType, libcommute::fermion> & MOp,
                         const StatesClassification &S,
                         const HamiltonianPart &HFrom,
                         const HamiltonianPart &HTo) :
      Complex(std::is_same<ScalarType, ComplexType>::value || HFrom.isComplex() || HTo.isComplex()),
      MOp(&MOp), S(S), HFrom(HFrom), HTo(HTo)
    {}

    /** Compute all the matrix elements. Changes the Status of the object to Computed. */
    void compute();

    void setFromAdjoint(const MonomialOperatorPart &part);

    /** Print all matrix elements of the operator to screen. */
    void print_to_screen() const;

    bool isComplex() const { return Complex; }

    /** Returns the row ordered sparse matrix of matrix elements. */
    template<bool Complex> RowMajorMatrixType<Complex>& getRowMajorValue();
    template<bool Complex> const RowMajorMatrixType<Complex>& getRowMajorValue() const;
    /** Returns the column ordered sparse matrix of matrix elements. */
    template<bool Complex> ColMajorMatrixType<Complex>& getColMajorValue();
    template<bool Complex> const ColMajorMatrixType<Complex>& getColMajorValue() const;
    /** Returns the right hand side index. */
    BlockNumber getRightIndex() const;
    /** Returns the left hand side index. */
    BlockNumber getLeftIndex() const;

    #ifdef ENABLE_SAVE_PLAINTEXT
    /** Save the data to the file.
     * \param[in] path Path to the file.
     */
    bool savetxt(const boost::filesystem::path &path);
    #endif

private:

    template<bool C, bool HC> void computeImpl();
    template<bool C> void print_to_screenImpl() const;
};

} // end of namespace Pomerol
#endif // endif :: #ifdef __INCLUDE_FIELDOPERATORPART_H