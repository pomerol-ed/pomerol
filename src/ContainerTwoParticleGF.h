#ifndef __CONTAINER2PGF__
#define __CONTAINER2PGF__

#include "config.h"

class TwoParticleGFContainer;


class TwoParticleGFContainer
{ 
  struct IndexCombination;
  std::vector<IndexCombination*> NonTrivialCombinations;
};

struct TwoParticleGFContainer::IndexCombination
{
  unsigned short Indices[4];
};
#endif // endif :: #ifndef __CONTAINER2PGF__
