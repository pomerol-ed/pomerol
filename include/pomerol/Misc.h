/** \file Misc.h
**    \brief Declares very common type names and macros.
**
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef __INCLUDE_MISC_H
#define __INCLUDE_MISC_H

#include<pomerol/first_include.h>

#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>
#include<iomanip>
#include<type_traits>

#include <libcommute/loperator/state_vector.hpp>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include<Eigen/Core>
#include<Eigen/Sparse>

#ifdef ENABLE_SAVE_PLAINTEXT
#include <boost/filesystem.hpp>
#endif

#ifdef POMEROL_USE_OPENMP
#include <omp.h>
#endif

namespace Pomerol {

#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#ifndef NDEBUG
#define DEBUG(MSG)            std::cout << MSG_PREFIX << MSG << std::endl
#else
#define DEBUG(MSG)            NULL;
#endif
#define INFO(MSG)             std::cout << MSG << std::endl
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush
#define ERROR(MSG)            std::cerr << MSG_PREFIX << MSG << std::endl

/** Real floating point type. */
typedef double RealType;
/** Complex type. */
typedef std::complex<double> ComplexType;

/** Index represents a combination of spin, orbital, and lattice indices **/
typedef unsigned int ParticleIndex;

// FIXME
//const FockState ERROR_FOCK_STATE = {}; // A state with the size==0 is an error state

/** Each Quantum State in the finite system is associated with a number.
 * This works for any basis, including Fock and Hamiltonian eigenbasis.
 * The Fock States are converted naturally from bitsets to ints.
 **/
using QuantumState = libcommute::sv_index_type;

/** Index represents a combination of spin, orbital, and lattice indices **/
typedef unsigned int ParticleIndex;

enum OperatorStatistics {fermion, boson};

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> LowerTriangularRealMatrixType;

template<bool Complex>
using MelemType = typename std::conditional<Complex, ComplexType, RealType>::type;

template<bool Complex>
using MatrixType = Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor>;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;
template<bool Complex>
using VectorType = Eigen::Matrix<MelemType<Complex>, Eigen::Dynamic,1,Eigen::AutoAlign>;

/** Sparse complex matrix */
template<bool Complex>
using ColMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::ColMajor>;
template<bool Complex>
using RowMajorMatrixType = Eigen::SparseMatrix<MelemType<Complex>, Eigen::RowMajor>;
//typedef Eigen::Triplet<RealType> RealTypeTriplet;
//typedef Eigen::Triplet<ComplexType> ComplexTypeTriplet;

/** Possible spin projections are \b down and \b up */
enum spin : unsigned short {down, up};

std::ostream & operator<<(std::ostream & os, spin s);

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems

/** Generalized 'square' function. */
template<typename T> inline T sqr(T x) { return x*x; }

/** */
#define DUMP_FLOATING_POINT_NUMBERS    10

//@{
/** Do-It-Once environment from A. Rubtsov
**
** When you want a piece of code to run exactly once, just write:
** \verbatim
do_once
    ... your code goes here...
end_do_once
\endverbatim
**/
#define do_once { static bool done_once=false; if (!done_once) {done_once=true;
#define end_do_once }; };
//@}

#define CHECK_MATSUBARA_NUM(num,num_of_matsubaras)	(num) < (num_of_matsubaras) && (num) >= -(num_of_matsubaras)

/** Permutation of 3 elements */
struct Permutation3 {
    const size_t perm[3];
    const int sign;
    bool operator==(const Permutation3& rhs) const;
    bool operator!=(const Permutation3& rhs) const;
    friend std::ostream& operator<<(std::ostream& out, const Permutation3 &rhs);
};
extern const Permutation3 permutations3[6];

/** Permutation of 4 elements */
struct Permutation4 {
    const size_t perm[4];
    const int sign;
    bool operator==(const Permutation4& rhs) const;
    bool operator!=(const Permutation4& rhs) const;
    friend std::ostream& operator<<(std::ostream& out, const Permutation4 &p);
};
extern const Permutation4 permutations4[24];

/** A tool to wrap the input and output of values. */
template <typename T> struct __num_format;
template <typename T> std::ostream& operator<<(std::ostream& lhs, const __num_format<T> &in);
template <typename T> std::istream& operator>>(std::istream& lhs, __num_format<T> &out);
template <typename T>
struct __num_format {
    static const int _prec = 12;
    T _v;
    __num_format(T v):_v(v){};
    operator T(){return _v;};
/*    void savetxt(const std::string& filename) {
        std::cout << "Saving " << typeid(*this).name() << " to " << filename << std::endl;
        std::ofstream out; out.open(filename.c_str()); out << *this << std::endl; out.close();
    };
*/
    friend std::ostream& operator<< <>(std::ostream& lhs, const __num_format<T> &in);
    friend std::istream& operator>> <>(std::istream& lhs, __num_format<T> &out);
};

template <typename T>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format<T> &in) {lhs << std::setprecision(in._prec) << in._v; return lhs;};
template <typename T>
inline std::istream& operator>>(std::istream& lhs, __num_format<T> &out) {lhs >> out._v; return lhs;};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format<ComplexType> &in){lhs << std::setprecision(in._prec) << real(in._v) << " " << imag(in._v); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, __num_format<ComplexType> &out){RealType re,im; lhs >> re; lhs >> im; out._v = re+I*im; return lhs;};

/*template <class T>
inline bool __is_zero(const T& in, RealType threshold = std::numeric_limits<RealType>::epsilon()){return (std::abs(in)<threshold);};
*/

} // end of namespace Pomerol

/**
 * \mainpage
 * The source code and fetch instructions are located at <a href="http://pomerol.googlecode.com">project's Google code page</a>.
 * \section   ref_API libpomerol API
 * The general sequence of a calculation is:
 * -    Define all indices of each mode of the system, i.e. site + spin indices ( by IndexClassification ).
 *      This is also a moment to discover symmetries against permutations of indices.
 * -    Define the Fock space - create QuantumStates and sort them into the blocks by their quantum numbers ( by StatesClassification ).
 * -    Enter blocks of the Hamiltonian ( by Hamiltonian and HamiltonianPart) and diagonalize them.
 * -    Find creation and annihilation operators in the eigenbasis of the Hamiltonian ( by FieldOperator and FieldOperatorPart ).
 * -    Calculate thermal quantities:
 *      -   The density matrix ( by DensityMatrix and DensityMatrixPart ).
 *      -   The Green's function ( by GreensFunction, GreensFunctionPart and GFContainer to store the values of GF for various index combinations ).
 *      -   The TwoParticle Greens Function ( by TwoParticleGF, TwoParticleGFPart and TwoParticleGFContainer ).
 * -    Calculcate the Vertex Function out of GF- and TwoParticleGF- containers - no work with Fock space is done ( by Vertex4 ).
 *
 * A hint: Refer to <a href="inherits.html">a Class Hierarchy</a> if provided
 *
 * \section ref_conventions Conventions
 *
 * \par Green's function
 * \f[
 *      G(\omega_n) = -\int_0^\beta \langle\mathbf{T}c_i(\tau)c^+_j(0)\rangle e^{i\omega_n\tau} d\tau
 * \f]
 *
 * \par Two-particle Green's function:
 * \f[ \chi_{ijkl}(\omega_{n_1},\omega_{n_2};\omega_{n_3},\omega_{n_1}+\omega_{n_2}-\omega_{n_3}) =
 *   \int_0^\beta
 *     \langle\mathbf{T} c_i(\tau_1)c_j(\tau_2)c^+_k(\tau_3)c^+_l(0) \rangle
 *     \exp(i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3)
 *   d\tau_1 d\tau_2 d\tau_3
 * \f]
 *
 * \par The Wick part of a two-particle Green's function:
 * \f[
 * \chi^0_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) =
 * \beta\delta_{\omega_1\omega_4}\delta_{\omega_2\omega_3}G_{14}(\omega_1)G_{23}(\omega_2) -
 *  \beta\delta_{\omega_1\omega_3}\delta_{\omega_2\omega_4}G_{13}(\omega_1)G_{24}(\omega_2)
 * \f]
 * \par An irreducible vertex part:
 * \f[ \Gamma_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) \equiv
 *     \chi_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) -
 *     \chi^{0}_{1234}(\omega_1,\omega_2;\omega_3,\omega_4)
 * \f]
 *
 * \par An amputated irreducible vertex part:
 * \f[ \gamma_{1234}(\omega_1,\omega_2;\omega_3,\omega_4) \equiv
 *     \sum_{1'2'3'4'}
 *     (G^{-1}(\omega_1))_{11'} (G^{-1}(\omega_2))_{22'}
 *     \Gamma_{1'2'3'4'}(\omega_1,\omega_2;\omega_3,\omega_4)
 *     (G^{-1}(\omega_3))_{3'3} (G^{-1}(\omega_4))_{4'4}
 * \f]
 *
 * \section Features
 * -    OpenMP support
 * \todo HDF5 storage ( read / write )
 * \section ref_authors Authors
 * -    Andrey Antipov <antipov[at]ct-qmc.org>
 * -    Igor Krivenko <igor[at]shg.ru>
 *
 * with a help from Mikhail Alejnikov, Alexey Rubtsov, Christoph Jung and Aljoscha Wilhelm
 *
 * We acknowledge <em>RRC Kurchatov Institute</em> for providing computing resources.
 *
 */

#endif // #ifndef __INCLUDE_MISC_H

