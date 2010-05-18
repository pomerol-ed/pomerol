// pomerol/trunk/src/Term.h
// This file is a part of pomerol diagonalization code

/** \file Term.h
**  \brief Declaration of Term<>, nnTerm, spinflipTerm and other term classes.
** 
**  \author	Andrey Antipov (antipov@ct-qmc.org)
*/


#ifndef __TERM_H__
#define __TERM_H__

#include "config.h"
#include <map>
#include <string>
#include <ostream>

// Forward declarations
template <int len> class Term;
template <int len> std::ostream& operator<<(std::ostream& output, const Term<len>& out);

/** Term<len> - a class to represent a term in a formula. 
 * The term is meant to be a finite number (defined by template parameter len) of creation and annihilation operators. 
 * Their indices and order are stored in this class. Examples are nn or spinflip or other type of terms.
 */
template <int len> class Term 
{
public:
   string type;
   bool diag; 			//!< Determine, whether the term will produce diagonal or non-diagonal matrix element.
   bool order[len]; 		//!< Order is a sequence of true/false values. True means operation operator, false - annihilation operator.
   unsigned short bit[len]; 	//!< The bits which are used by field operators.
   RealType Value;		//!< Define which value should be added to Hamiltonian by this term. Actually, this is the matrix element.
   friend std::ostream & (operator << <>) (ostream &, const Term<len> &); //!< Output of the class 
};

/** nnTerm - class, inherited from Term<4> class. Represents an nn type of term.
 * This type of term is diagonal and has two pairs of equal indices. The order is creation-annihilation-creation-annihilation operators.
 */
class nnTerm : public Term<4> {
public:
   nnTerm(unsigned short bit1, unsigned short bit2, RealType Val);
   friend std::ostream & (operator <<) (ostream &,const nnTerm&);
};

/** spinflipTerm - class, inherited from Term<4> class. Represents a spinflip type of term.
 * This type of term is non-diagonal and has the defined order of two creation operators followed by two annihilation operators.
 */
class spinflipTerm : public Term<4> {
public:
   spinflipTerm(unsigned short bit1, unsigned short bit2, unsigned short bit3, unsigned short bit4, RealType Val);
   friend std::ostream & (operator <<) (ostream &,const spinflipTerm&);
};

/*==A GENERIC TERM METHODS======================================*/

template <int len> std::ostream& operator<<(std::ostream& output,const Term<len>& out)
{
if (len == 4 && out.type == "nn") { output << (nnTerm&) out; return output; }; 

     output << "A ";
     output << ((out.diag)?"diagonal ":"non-diagonal ");
     output << out.type << " term : " << out.Value << "*(";

     for (unsigned short it=0;it<len;++it)
     {
        output<<"c" << (out.order[it]?"^+":"") << "_{" << out.bit[it] << "}";
     };

     output << ")";
     return output;
}
/*
template <int len> Term<len>::~Term() {
	delete[] order;
	delete[] bit;
};
*/
#endif
