/** \file Misc.h
**    \brief Declares very common type names and macros.
** 
** \author    Igor Krivenko (igor@shg.ru)
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

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include<Eigen/Core>
#include<Eigen/Sparse>

/** Real floating point type. */
typedef double RealType;
/** Complex type. */
typedef std::complex<RealType> ComplexType;

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign> MatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign> RealMatrixType;

typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign> LowerTriangularRealMatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> VectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;

/** Sparse complex matrix */
typedef Eigen::SparseMatrix<RealType,Eigen::ColMajor> ColMajorMatrixType;
typedef Eigen::SparseMatrix<RealType,Eigen::RowMajor> RowMajorMatrixType;
typedef Eigen::DynamicSparseMatrix<RealType,Eigen::ColMajor> DynamicSparseMatrixType;

/** Index represents a combination of spin, orbital, and lattice indices **/
typedef unsigned short ParticleIndex;

/** Possible spin projections are \b down and \b up */
enum spin {down, up};

/** This represent the progress of calculation of any complex object in the code */
enum ObjectStatus {Constructed, Prepared, Computed};

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems


/** Generalized 'square' function. */
template<typename T> inline T sqr(T x) { return x*x; }

/** Easy enumeration for orbital names */
enum OrbitalValue {s=0, p=1, d=2, f=3};             //!< The enum for s,p,d,f - orbitals

/** Permutation of 3 elements */
struct Permutation3 {
    const size_t perm[3];
    const int sign;
};

/** Permutation of 4 elements */
struct Permutation4 {
    const size_t perm[4];
    const int sign;
};

static const Permutation4 permutations4[24] = {
    {{0,1,2,3}, 1},  {{0,1,3,2},-1},  {{0,2,1,3},-1},  {{0,2,3,1}, 1},  {{0,3,1,2}, 1},  {{0,3,2,1},-1},
    {{1,0,2,3},-1},  {{1,0,3,2}, 1},  {{1,2,0,3}, 1},  {{1,2,3,0},-1},  {{1,3,0,2},-1},  {{1,3,2,0}, 1},
    {{2,0,1,3}, 1},  {{2,0,3,1},-1},  {{2,1,0,3},-1},  {{2,1,3,0}, 1},  {{2,3,0,1}, 1},  {{2,3,1,0},-1},
    {{3,0,1,2},-1},  {{3,0,2,1}, 1},  {{3,1,0,2}, 1},  {{3,1,2,0},-1},  {{3,2,0,1},-1},  {{3,2,1,0}, 1}
};
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

#define INFO(MSG)                 std::cout << MSG << std::endl;
#define INFO_NONEWLINE(MSG)       std::cout << MSG << std::flush;
#define ERROR(MSG)                std::cerr << MSG << std::endl;

#ifdef NDEBUG
  #define DEBUG(x)
#else
  #define DEBUG(x)  INFO(x)
#endif

#define num_cout std::cout << std::setprecision(8) << std::setw(9) << std::left
#define iomanip_prefs std::setprecision(8) << std::setw(9) << std::left

#endif // #ifndef __INCLUDE_MISC_H
