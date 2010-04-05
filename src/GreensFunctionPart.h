#ifndef ____DEFINE_GREENS_FUNCTION_PART____
#define ____DEFINE_GREENS_FUNCTION_PART____
#include "config.h"
#include "getStates.h"
#include "hpart.h"
#include "CCXpart.h"
#include "FieldOperators.h"
#include "output.h"

class GreensFunctionPart
{
  int w;
  int i;
  int j;
  getStates &System;
  FieldOperatorPart &CX;
  FieldOperatorPart &C;
  getHpart &H_n;
  getHpart &H_m;
  output_handle OUT;
  
  RealType value;

  public:
    GreensFunctionPart(int w_, int i_, int j_, getStates &S_, FieldOperatorPart &CX_, FieldOperatorPart &C_, getHpart &H_n_, getHpart &H_m_, output_handle &OUT_):
	    w(w_),i(i_),j(j_),System(S_), CX(CX_), C(C_), H_n(H_n_), H_m(H_m_), OUT(output_handle(OUT_.path()+"//Gw")){};
    RealType& getValue();
};

#endif
