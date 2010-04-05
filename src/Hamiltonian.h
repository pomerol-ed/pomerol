#ifndef ____DEFINE_HAMILTONIAN____
#define ____DEFINE_HAMILTONIAN____
#include "config.h"
#include "BitClassification.h"
#include "getStates.h"
#include "hpart.h"
#include "output.h"
#include <vector>

class Hamiltonian
{
  getHpart **Hpart;


  BitClassification &Formula;
  getStates& S;
  output_handle &OUT;
  string config_path;

public :

  Hamiltonian(BitClassification &F_, getStates &S_,output_handle &OUT_, string &config_path_);
  void enter(bool diag=false,bool dump=false);

  getHpart& block(const QuantumNumbers &in);
  getHpart& block(BlockNumber in);
  RealType eigenval( QuantumState &state );

  void diagonalize();
  void dump();
};

#endif // endif :: #ifndef ____DEFINE_HAMILTONIAN____

