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

/** \file Misc.h
**    \brief Declares very common type names and macros.
** 
** \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author    Alexey Rubtsov (alex@shg.ru)
** \author    Andrey Antipov (antipov@shg.ru)
*/

#ifndef __INCLUDE_MISC_H
#define __INCLUDE_MISC_H

/** \cond */
//#define NDEBUG
/** \endcond */


#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>
#include<iomanip>


#include<boost/shared_ptr.hpp>
#include<boost/scoped_ptr.hpp>
#include<boost/make_shared.hpp>
#include<boost/dynamic_bitset.hpp>
#include<boost/tuple/tuple.hpp>
#include<boost/utility.hpp>


#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include<Eigen/Core>
#include<Eigen/Sparse>

#ifdef ENABLE_SAVE_PLAINTEXT
#include <boost/filesystem.hpp>
#endif

#include <boost/mpi.hpp>

#define REALTYPE_DOUBLE

namespace Pomerol{

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

/** Matrix element type. */
#ifdef POMEROL_COMPLEX_MATRIX_ELEMENS
#warning using complex matrix elements
typedef ComplexType MelemType;
#else
typedef RealType MelemType;
#endif

/** Index represents a combination of spin, orbital, and lattice indices **/
typedef unsigned int ParticleIndex;

//enum Statistics { fermion, boson };
//typedef boost::tuple<bool,Statistics,ParticleIndex> AtomicOp;
//typedef AtomicOp<1,fermion,ParticleIndex> AtomicCdag;

/** Fock State representation. */ 
typedef boost::dynamic_bitset<> FockState;
const FockState ERROR_FOCK_STATE = FockState(); // A state with the size==0 is an error state

/** Each Quantum State in the finite system is associated with a number. 
 * This works for any basis, including Fock and Hamiltonian eigenbasis. 
 * The Fock States are converted naturally from bitsets to ints. 
 **/
typedef unsigned long QuantumState;

/** Index represents a combination of spin, orbital, and lattice indices **/
typedef unsigned int ParticleIndex;

enum OperatorStatistics {fermion, boson}; 
/** A creation and annihilation operators */
//typedef boost::tuple<bool,OperatorStatistics,ModeIndex,ParticleIndex> ElementaryOperator; 
//typedef std::tuple<bool,OperatorStatistics,ParticleIndex> ElementaryOperator; 
//template<ParticleIndex P> using ElemCreatOpFermion = typename ElementaryOperator<1,0,P>;
//template<ParticleIndex> typedef boost::tuple<0,0,ParticleIndex> ElemAnnihOpFermion; 

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> LowerTriangularRealMatrixType;
/** Default Matrix Type comes from MelemType. */
typedef Eigen::Matrix<MelemType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> MatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;
/** Default vector type comes from MelemType. */
typedef Eigen::Matrix<MelemType,Eigen::Dynamic,1,Eigen::AutoAlign> VectorType;

/** Sparse complex matrix */
typedef Eigen::SparseMatrix<MelemType,Eigen::ColMajor> ColMajorMatrixType;
typedef Eigen::SparseMatrix<MelemType,Eigen::RowMajor> RowMajorMatrixType;
typedef Eigen::DynamicSparseMatrix<MelemType,Eigen::ColMajor> DynamicSparseMatrixType;
//typedef Eigen::Triplet<RealType> RealTypeTriplet;
//typedef Eigen::Triplet<ComplexType> ComplexTypeTriplet;

/** Possible spin projections are \b down and \b up */
enum spin {down, up};

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
 * \section ref_exec pomerolDiag       
 * \verbatim pomerolDiag -b <beta> -m <Nm> -l <LatticeFile> \endverbatim
 * \param beta The value of inverse temperature for Thermal objects (i.e DensityMatrix, TwoParticleGF, GreensFunction etc).     
 * \param Nm A number of positive Matsubara frequencies, for which the computation of TwoParticleGF and Vertex4 should be performed. The GreensFunction is obtained with 8*Nm Matsubara Frequencies.      
 * \param LatticeFile A lattice .json file, see page \ref ref_lattice .
 *
 * Call "pomerolDiag -h" for this information
 * \section   ref_API libpomerol API
 * The general sequence of a calculation is:
 * -    Open a Lattice .json file ( by LatticeAnalysis ).
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
 *   -\int_0^\beta
 *      \langle\mathbf{T} c_i(\tau_1)c_j(\tau_2)c^+_k(\tau_3)c^+_l(0) \rangle
 *      \exp(i\omega_{n_1}\tau_1+i\omega_{n_2}\tau_2-i\omega_{n_3}\tau_3)
 *    d\tau_1 d\tau_2 d\tau_3
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

