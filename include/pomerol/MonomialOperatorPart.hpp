/** \file include/pomerol/FieldOperatorPart.h
** \brief Declaration of FieldOperatorPart, CreationOperatorPart and AnnihilationOperatorPart classes.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_MONOMIALOPERATORPART_H
#define POMEROL_INCLUDE_MONOMIALOPERATORPART_H

#include "HamiltonianPart.hpp"
#include "HilbertSpace.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <memory>
#include <ostream>
#include <type_traits>

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

private:
    bool MOpComplex;
    void const* MOp;

    bool Complex;

    /** A reference to the StateClassification object. */
    StatesClassification const& S;
    /** A reference to the HamiltonianPart on the right hand side. */
    HamiltonianPart const& HFrom;
    /** A reference to the HamiltonianPart on the left hand side. */
    HamiltonianPart const& HTo;

protected:
    /** Storage of the matrix elements of the operator. Row ordered sparse matrix. */
    std::shared_ptr<void> elementsRowMajor;
    /** Copy of the Storage of the matrix elements of the operator. Column ordered sparse matrix. */
    std::shared_ptr<void> elementsColMajor;
    /** The tolerance with which the matrix elements are evaluated. */
    RealType const MatrixElementTolerance = 1e-8;

public:
    /** Constructor.
     * \param[in] IndexInfo A const reference to the IndexClassification object.
     * \param[in] S A const reference to the StateClassification object.
     * \param[in] HFrom A const reference to the HamiltonianPart on the right hand side.
     * \param[in] HTo A const reference to the HamiltonianPart on the left hand side.
     * \param[in] PIndex Index of the field operator.
     */
    template <typename ScalarType>
    MonomialOperatorPart(LOperatorType<ScalarType> const& MOp,
                         StatesClassification const& S,
                         HamiltonianPart const& HFrom,
                         HamiltonianPart const& HTo)
        : MOpComplex(std::is_same<ScalarType, ComplexType>::value),
          MOp(&MOp),
          Complex(MOpComplex || HFrom.isComplex() || HTo.isComplex()),
          S(S),
          HFrom(HFrom),
          HTo(HTo) {}

    /** Compute all the matrix elements. Changes the Status of the object to Computed. */
    void compute();

    void setFromAdjoint(MonomialOperatorPart const& part);

    bool isComplex() const { return Complex; }

    /** Returns the row ordered sparse matrix of matrix elements. */
    template <bool C> RowMajorMatrixType<C>& getRowMajorValue();
    template <bool C> RowMajorMatrixType<C> const& getRowMajorValue() const;
    /** Returns the column ordered sparse matrix of matrix elements. */
    template <bool C> ColMajorMatrixType<C>& getColMajorValue();
    template <bool C> ColMajorMatrixType<C> const& getColMajorValue() const;
    /** Returns the right hand side index. */
    BlockNumber getRightIndex() const { return HFrom.getBlockNumber(); }
    /** Returns the left hand side index. */
    BlockNumber getLeftIndex() const { return HTo.getBlockNumber(); }

    friend std::ostream& operator<<(std::ostream& os, MonomialOperatorPart const& part) {
        if(part.isComplex())
            part.streamOutputImpl<true>(os);
        else
            part.streamOutputImpl<false>(os);
        return os;
    }

private:
    template <bool C, bool HC> void computeImpl();
    template <bool C> void streamOutputImpl(std::ostream& os) const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_MONOMIALOPERATORPART_H
