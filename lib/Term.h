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


// pomerol/trunk/src/Term.h
// This file is a part of pomerol diagonalization code

/** \file Term.h
**  \brief Declaration of Term<>, nnTerm, spinflipTerm and other term classes.
** The manifold of terms describes a formula which defines Hamiltonian on a given lattice
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/


#ifndef __INCLUDE_TERM_H
#define __INCLUDE_TERM_H

#include"Misc.h"

/** Term - a class to represent a term in a formula. 
 * The term is meant to be a finite number (defined by template parameter len) of creation and annihilation operators. 
 * Their indices and order are stored in this class. Examples are nn or spinflip or other type of terms.
 */
class Term 
{
public:
   int N;            //!< Amount of field operators in term
   std::string type;
   bool diag;             //!< Determine, whether the term will produce diagonal or non-diagonal matrix element.
   bool *order;         //!< Order is a sequence of true/false values. True means creation operator, false - annihilation operator.
   unsigned short *bit;     //!< The bits which are used by field operators.
   RealType Value;        //!< Define which value should be added to Hamiltonian by this term. Actually, this is the matrix element.
   friend std::ostream & (operator <<) (std::ostream &, const Term &); //!< Output of the class 
   virtual ~Term();
};

/** nnTerm - class, inherited from Term class with N=4. Represents an nn type of term.
 * This type of term is diagonal and has two pairs of equal indices. The order is creation-annihilation-creation-annihilation operators.
 */
class nnTerm : public Term {
public:
   nnTerm(unsigned short bit1, unsigned short bit2, RealType Val);
   friend std::ostream & (operator <<) (std::ostream &,const nnTerm&);
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
#endif // #define __INCLUDE_TERM_H
