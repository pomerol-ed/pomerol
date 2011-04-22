#include "Term.h"

nnTerm::nnTerm(unsigned short bit1, unsigned short bit2, RealType Val)
{
  diag = true;
  Value = Val;
  N=4;
  bit = new unsigned short [4];
  bit[0]=bit1; bit[1]=bit1;
  bit[2]=bit2; bit[3]=bit2;
  
  order = new bool [4];
  order[0]=true;
  order[1]=false;
  order[2]=true;
  order[3]=false;

  type="nn";
}

std::ostream& operator<<(std::ostream& output, const nnTerm& out)
{
   output << "Diagonal " << out.type << " term, " << out.Value << "*(n_{" <<out.bit[0] <<"}n_{"<< out.bit[2] << "})"; 
   return output;
}

spinflipTerm::spinflipTerm(unsigned short bit1, unsigned short bit2, unsigned short bit3, unsigned short bit4, RealType Val)
{
  diag = false;
  Value = Val;
  N=4;
  bit = new unsigned short [4];
  bit[0]=bit1; bit[1]=bit2;
  bit[2]=bit3; bit[3]=bit4;
  
  order = new bool [4];
  order[0]=true;
  order[1]=true;
  order[2]=false;
  order[3]=false;

  type="spinflip";
}

nTerm::nTerm(unsigned short bit1, RealType Val)
{
  diag = true;
  Value = Val;
  N=2;
  bit = new unsigned short [2];
  bit[0]=bit1; bit[1]=bit1;
  
  order = new bool [2];
  order[0]=true;
  order[1]=false;

  type="n";
};

std::ostream& operator<<(std::ostream& output,const Term& out)
{
if (out.N == 4 && out.type == "nn") { output << (nnTerm&) out; return output; }; 

     output << "A ";
     output << ((out.diag)?"diagonal ":"non-diagonal ");
     output << out.type << " term : " << out.Value << "*(";

     for (unsigned short it=0;it<out.N;++it)
     {
        output<<"c" << (out.order[it]?"^+":"") << "_{" << out.bit[it] << "}";
     };

     output << ")";
     return output;
}

Term::~Term() {
	delete[] order;
	delete[] bit;
};


