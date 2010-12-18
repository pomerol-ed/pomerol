#ifndef ____VERTEX4____
#define ____VERTEX4____

#include "GreensFunction.h"
#include "TwoParticleGF.h"

class Vertex4 {

    TwoParticleGF &Chi;
    GreensFunction &g1;
    GreensFunction &g2;
    GreensFunction &g3;
    GreensFunction &g4;
    
    unsigned short ChiBit1, ChiBit2, ChiBit3, ChiBit4;
  
public:
    
    Vertex4(TwoParticleGF &Chi, GreensFunction &g1, GreensFunction &g2, GreensFunction &g3, GreensFunction &g4);
    
    ComplexType operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);
    ComplexType getAmputated(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3);

    //void dumpMatsubara(unsigned short points);
    //void dumpAmputatedMatsubara(unsigned short points);
};

#endif // endif :: #ifndef ____VERTEX4____