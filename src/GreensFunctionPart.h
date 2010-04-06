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
    GreensFunctionPart(getStates &S_, AnnihilationOperatorPart& C, CreationOperatorPart& CX, output_handle &OUT):
	    System(S), CX(CX), C(C), OUT(output_handle(OUT.path()+"//Gw")){};
    RealType& getValue();
};

#endif
