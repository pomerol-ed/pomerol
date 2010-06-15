#ifndef ____DEFINE_HAMILTONIAN____
#define ____DEFINE_HAMILTONIAN____
#include "config.h"
#include "BitClassification.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "output.h"
#include <vector>

class Hamiltonian
{
  HamiltonianPart **Hpart;


  BitClassification &Formula;
  StatesClassification& S;
  output_handle &OUT;
  string config_path;

  RealType GroundEnergy;
public :

  Hamiltonian(BitClassification &F_, StatesClassification &S_,output_handle &OUT_, string &config_path_);
  void enter();

  HamiltonianPart& part(const QuantumNumbers &in);
  HamiltonianPart& part(BlockNumber in);
  RealType eigenval( QuantumState &state );
  RealType getGroundEnergy();

  void diagonalize();
  void dump();
  void reduce(const RealType Cutoff);
private:
  void computeGroundEnergy();
};

#endif // endif :: #ifndef ____DEFINE_HAMILTONIAN____

