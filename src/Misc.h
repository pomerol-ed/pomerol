//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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

#ifdef pomerolOpenMP
#include <omp.h>
#endif

#define REALTYPE_DOUBLE

namespace Pomerol{

/** Real floating point type. */
typedef double RealType;
/** Complex type. */
typedef std::complex<RealType> ComplexType;

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> MatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> LowerTriangularRealMatrixType;

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

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems

/** Generalized 'square' function. */
template<typename T> inline T sqr(T x) { return x*x; }

/** Easy enumeration for orbital names. */
enum OrbitalValue {s=0, p=1, d=2, f=3};             //!< The enum for s,p,d,f - orbitals

/** Permutation of 3 elements */
struct Permutation3 {
    const size_t perm[3];
    const int sign;
    bool operator==(const Permutation3& rhs) const {return (sign==rhs.sign && perm[0] == rhs.perm[0] && perm[1]==rhs.perm[1]);};
    bool operator!=(const Permutation3& rhs) const {return !(*this==rhs);};
    friend std::ostream& operator<<(std::ostream& out, const Permutation3 &rhs)
        { out << (rhs.sign==-1?"-":" ") << rhs.perm[0]+1 << rhs.perm[1]+1 << rhs.perm[2]+1 << std::flush; return out;}; 
};

/** Permutation of 4 elements */
struct Permutation4 {
    const size_t perm[4];
    const int sign;
    bool operator==(const Permutation4& rhs) const {return (sign==rhs.sign && perm[0] == rhs.perm[0] && perm[1]==rhs.perm[1] && perm[2] == rhs.perm[2]);};
    bool operator!=(const Permutation4& rhs) const {return !(*this==rhs);};
    friend std::ostream& operator<<(std::ostream& out, const Permutation4 &rhs)
        { out << (rhs.sign==-1?"-":" ") << rhs.perm[0]+1 << rhs.perm[1]+1 << rhs.perm[2]+1 << rhs.perm[3]+1 << std::flush; return out;}; 
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

#define CHECK_MATSUBARA_NUM(num,num_of_matsubaras)	(num) < (num_of_matsubaras) && (num) >= -(num_of_matsubaras)

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
 * \page ref_lattice Lattice JSON-file syntax
 * \section ref_json_generic Generic Information
 * Several examples for json Lattice files are provided in \b doc/lattices subdirectory. A python script, generating such 
 * file is also provided in \b scripts subdirectory: \a VLALatticeGenerator.py .
 *
 * For a generic information about JSON config files, refer to a <a href="http://en.wikipedia.org/wiki/JSON"> wiki page</a>.
 *
 * \section ref_json_examples Examples
 * \subsection ref_json_2site An example of a Hubbard model on a 2-site cluster 
 *
 * \warning at this point no other models are supported.
 *
 * A lattice file should contain information about each site of the lattice, they are named by their numbers
 * \param type Orbital complexity of current site. A support of "p"-orbital sites is experimental and subject to strong testing. For this reason no multiorbital example is provided.
 * \param U A coulomb repulsion between electrons on a site.
 * \param LocalMu A chemical potential.
 * \param hopping A list of hoppings to the other sites.
 * \par Hopping parameters
 * \param "to" A site which is connected by a hopping integral.
 * \param "value" Value of hopping integral.
 * \param "orbital_from" Orbital index at local site, where hopping occurs from. For s-orbital case this parameter is always 0.
 * \param "orbital_to" Orbital index at remote site, where hopping occurs to. For s-orbital case this parameter is always 0.
 * \verbinclude lattices/Lattice2.json
 */

#endif // #ifndef __INCLUDE_MISC_H

