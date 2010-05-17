#ifndef __TERM_H__
#define __TERM_H__

#include "config.h"
#include <map>
#include <string>
#include <ostream>

/*== A GENERIC TERM=============================================*/

template <int len> class Term;
template <int len> std::ostream& operator<<(std::ostream& output, const Term<len>& out);

template <int len> class Term 
/* A class to represent a term in a formula
 * It can be either nn term, or spin_flip or other type of terms
 */
{
public:
   string type;
   bool diag; // Determine, whether the term is diagonal or not
   bool order[len]; // Order is a sequence of true/false values. True means operation operator, false - annihilation operator
   unsigned short bit[len]; // The bits which are
   RealType Value; // Define which value to add
  // ~Term();
   friend std::ostream & (operator << <>) (ostream &, const Term<len> &);  
};

/*==SPECIFIC TERMS=============================================*/

class nnTerm : public Term<4> {
public:
   nnTerm(unsigned short bit1, unsigned short bit2, RealType Val);
   friend std::ostream & (operator <<) (ostream &,const nnTerm&);
};

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
