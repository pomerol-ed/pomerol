/** \file include/pomerol/HamiltonianPart.h
** \brief Declaration of the HamiltonianPart class.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_HAMILTONIANPART_H
#define __INCLUDE_HAMILTONIANPART_H

#include "Misc.h"
#include "IndexClassification.h"
#include "StatesClassification.h"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <memory>
#include <type_traits>

#ifdef ENABLE_SAVE_PLAINTEXT
#include <boost/filesystem/path.hpp>
#endif

namespace Pomerol{

/** HamiltonianPart is a class, which stores and diagonalizes the block of the Hamiltonian, which corresponds to a set of given quantum numbers. */
class HamiltonianPart : public ComputableObject {

    bool Complex;

    template<bool C>
    using LOperatorType = libcommute::loperator<MelemType<C>, libcommute::fermion>;
    const void* HOp;

    const StatesClassification &S;

    /** The number of Block. Defined in StatesClassification. */
    BlockNumber Block;

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
    HamiltonianPart(const libcommute::loperator<ScalarType, libcommute::fermion> &HOp,
                    const StatesClassification &S,
                    BlockNumber Block) :
        Complex(std::is_same<ScalarType, ComplexType>::value),
        HOp(&HOp), S(S), Block(Block)
    {}

    /** Fill in the H matrix. */
    void prepare();
    /** Diagonalize the H matrix and get EigenValues. */
    void compute();

    bool reduce(RealType ActualCutoff); // Useless now

    bool isComplex() const { return Complex; }

    /** Return the total dimensionality of the H matrix. This corresponds to the one in StatesClassfication. */
    InnerQuantumState getSize() const;

    /** Get the eigenvalue of the H matrix.
     * \param[in] Number of eigenvalue. */
    RealType getEigenValue(InnerQuantumState state) const;

    /** Returns calculated eigenvalues. */
    const RealVectorType& getEigenValues() const;

    /** Return the hamiltonian part matrix. */
    template<bool Complex>
    const MatrixType<Complex>& getMatrix() const;
    template<bool Complex>
    MatrixType<Complex>& getMatrix();

    /** Return the lowest Eigenvalue of the current part. */
    RealType getMinimumEigenvalue() const;
    /** Return the eigenstate of the H matrix.
     * \param[in] Number of eigenvalue. */
    template<bool Complex>
    VectorType<Complex> getEigenState(InnerQuantumState state) const;

    /** Return the BlockNumber associated with the Hamiltonian part. */
    BlockNumber getBlockNumber() const { return Block; }

    /** Print the part of hamiltonian to screen. */
    void print_to_screen() const;

    /** Save the data to the file.
     * \param[in] path Path to the file.
     */
    #ifdef ENABLE_SAVE_PLAINTEXT
    bool savetxt(const boost::filesystem::path &path1);
    #endif

private:

    template<bool Complex> void initHMatrix();
    template<bool Complex> void prepareImpl();
    template<bool Complex> void computeImpl();
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIANPART_H
