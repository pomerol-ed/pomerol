//	werner-ng/src/config.h
//	This file is a part of 'Werner NG' project.
/** \file config.h
**	\brief Declares very common type names and macros.
** 
** \author	Igor Krivenko (igor@shg.ru)
** \author	Alexey Rubtsov (alex@shg.ru)
** \author	Andrey Antipov (antipov@shg.ru)
*/

#ifndef _CONFIG_H
#define _CONFIG_H

/** \cond */
//#define NDEBUG
/** \endcond */

#include<Eigen/Core>
#include<Eigen/Sparse>
#include<Eigen/QR>
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
typedef Eigen::SparseMatrix<RealType,Eigen::ColMajor> SparseMatrixType;
typedef Eigen::DynamicSparseMatrix<RealType,Eigen::ColMajor> DynamicSparseMatrixType;

/** Possible spin projections are \b down and \b up */
enum spin {down, up};

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);	// 'static' to prevent linking problems

/** Generalized 'square' function. */
template<typename T> inline T sqr(T x) { return x*x; }

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

/** In some functions matrix elements less then this value are treated as zero.*/
#define MATRIX_ELEMENT_TOLERANCE	1e-10
#define DUMP_FLOATING_POINT_NUMBERS	10

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

// Maybe temporary
using std::cout;
using std::endl;
using std::flush;
using std::complex;
using std::ostream;
using std::ofstream;
using std::vector;
using std::string;

#define INFO(MSG)       std::cout << MSG << std::endl;
#define ERROR(MSG)      std::cerr << MSG << std::endl;

#ifdef NDEBUG
  #define DEBUG(x)
#else
  #define DEBUG(x)  INFO(x)
#endif

#endif /* #ifndef _CONFIG_H */
