// pomerol/trunk/src/Term.h
// This file is a part of pomerol diagonalization code

/** \file Term.h
**  \brief Declaration of Term<>, nnTerm, spinflipTerm and other term classes.
** The manifold of terms describes a formula which defines Hamiltonian on a given lattice
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/


#ifndef __TERM_H__
#define __TERM_H__

#include "config.h"
#include <map>
#include <string>
#include <ostream>

/** Term - a class to represent a term in a formula. 
 * The term is meant to be a finite number (defined by template parameter len) of creation and annihilation operators. 
 * Their indices and order are stored in this class. Examples are nn or spinflip or other type of terms.
 */
class Term 
{
public:
   int N;            //!< Amount of field operators in term
   string type;
   bool diag;             //!< Determine, whether the term will produce diagonal or non-diagonal matrix element.
   bool *order;         //!< Order is a sequence of true/false values. True means creation operator, false - annihilation operator.
   unsigned short *bit;     //!< The bits which are used by field operators.
   RealType Value;        //!< Define which value should be added to Hamiltonian by this term. Actually, this is the matrix element.
   friend std::ostream & (operator <<) (ostream &, const Term &); //!< Output of the class 
   virtual ~Term();
};

/** nnTerm - class, inherited from Term class with N=4. Represents an nn type of term.
 * This type of term is diagonal and has two pairs of equal indices. The order is creation-annihilation-creation-annihilation operators.
 */
class nnTerm : public Term {
public:
   nnTerm(unsigned short bit1, unsigned short bit2, RealType Val);
   friend std::ostream & (operator <<) (ostream &,const nnTerm&);
};

/** spinflipTerm - class, inherited from Term class with N=4. Represents a spinflip type of term.
 * This type of term is non-diagonal and has the defined order of two creation operators followed by two annihilation operators.
 */
class spinflipTerm : public Term {
public:
   spinflipTerm(unsigned short bit1, unsigned short bit2, unsigned short bit3, unsigned short bit4, RealType Val);
};

class nTerm : public Term {
public:
   nTerm(unsigned short bit,RealType Val);
};
#endif
