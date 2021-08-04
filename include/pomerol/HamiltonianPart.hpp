/** \file include/pomerol/HamiltonianPart.h
** \brief Declaration of the HamiltonianPart class.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/
#ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_H
#define POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_H

#include "ComputableObject.hpp"
#include "IndexClassification.hpp"
#include "Misc.hpp"
#include "StatesClassification.hpp"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <iostream>
#include <memory>
#include <type_traits>

namespace Pomerol {

/** HamiltonianPart is a class, which stores and diagonalizes the block of the Hamiltonian, which corresponds to a set of given quantum numbers. */
class HamiltonianPart : public ComputableObject {

    StatesClassification const& S;

    /** The number of Block. Defined in StatesClassification. */
    BlockNumber Block;

    bool Complex;

    void const * HOp;

    /** A matrix filled with matrix elements of HamiltonianPart in the space of FockState's.
     *  After diagonalization it stores the eigenfunctions of the problem in a rows of H. */
    std::shared_ptr<void> HMatrix = nullptr;

    /** A vector of eigenvalues of the HamiltonianPart. */
    RealVectorType Eigenvalues;

    friend class Hamiltonian;

public:

    /** Constructor.
     * \param[in] IndexInfo IndexClassification object. Provides information about the indices in the problem.
     * \param[in] F IndexHamiltonian object. Provides all Terms required to fill the Hamiltonian part.
     * \param[in] S StatesClassification object. Provides information about Fock States of the problem.
     * \param[in] Block The BlockNumber of current part. It is a genuine id of the part. */
    template<typename ScalarType>
    HamiltonianPart(LOperatorType<ScalarType> const& HOp,
                    StatesClassification const& S,
                    BlockNumber Block) :
         S(S), Block(Block), Complex(std::is_same<ScalarType, ComplexType>::value), HOp(&HOp)
    {}

    /** Fill in the H matrix. */
    void prepare();
    /** Diagonalize the H matrix and get EigenValues. */
    void compute();

    bool reduce(RealType ActualCutoff);

    bool isComplex() const { return Complex; }

    /** Return the total dimensionality of the H matrix. This corresponds to the one in StatesClassfication. */
    InnerQuantumState getSize() const;

    /** Get the eigenvalue of the H matrix.
     * \param[in] Number of eigenvalue. */
    RealType getEigenValue(InnerQuantumState state) const;

    /** Returns calculated eigenvalues. */
    RealVectorType const& getEigenValues() const;

    /** Return the hamiltonian part matrix. */
    template<bool Complex> MatrixType<Complex> const& getMatrix() const;
    template<bool Complex> MatrixType<Complex>& getMatrix();

    /** Return the lowest Eigenvalue of the current part. */
    RealType getMinimumEigenvalue() const;

    /** Return the eigenstate of the H matrix.
     * \param[in] Number of eigenvalue. */
    template<bool Complex>
    VectorType<Complex> getEigenState(InnerQuantumState state) const;

    /** Return the BlockNumber associated with the Hamiltonian part. */
    BlockNumber getBlockNumber() const { return Block; }

    /** Print the part of hamiltonian to screen. */
    friend std::ostream & operator<<(std::ostream & os, HamiltonianPart const& part)
    {
        if(part.isComplex())
            os << part.getMatrix<true>() << std::endl;
        else
            os << part.getMatrix<false>() << std::endl;

        return os;
    }

private:

    template<bool C> void initHMatrix();
    template<bool C> void prepareImpl();
    template<bool C> void computeImpl();

    void checkComputed() const;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_POMEROL_HAMILTONIANPART_H
