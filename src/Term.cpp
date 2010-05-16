#include "Term.h"

nnTerm::nnTerm(unsigned short bit1, unsigned short bit2, RealType Val)
{
  diag = true;
  Value = Val;
  bit[0]=bit1; bit[1]=bit1;
  bit[2]=bit2; bit[3]=bit2;
  
  order[0]=true;
  order[1]=false;
  order[2]=true;
  order[3]=false;

  type="nn";
}

std::ostream& operator<<(std::ostream& output,const nnTerm& out)
{
   output << "Diagonal " << out.type << " term, " << out.Value << "*(n_{" <<out.bit[0] <<"}n_{"<< out.bit[2] << "})"; 
   return output;
}
