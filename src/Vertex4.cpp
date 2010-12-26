#include "Vertex4.h"

Vertex4::Vertex4(TwoParticleGF &Chi, GreensFunction &g1, GreensFunction &g2, GreensFunction &g3, GreensFunction &g4) :
Chi(Chi), g1(g1), g2(g2), g3(g3), g4(g4)
{
    ChiBit1 = Chi.getBit(0);
    ChiBit2 = Chi.getBit(1);
    ChiBit3 = Chi.getBit(2);
    ChiBit4 = Chi.getBit(3);

    if((ChiBit1 != g1.getBit(0) || ChiBit1 != g1.getBit(1)) ||
       (ChiBit2 != g2.getBit(0) || ChiBit1 != g2.getBit(1)) ||
       (ChiBit3 != g3.getBit(0) || ChiBit1 != g3.getBit(1)) ||
       (ChiBit4 != g4.getBit(0) || ChiBit1 != g4.getBit(1))
      ) assert(0);
}

ComplexType Vertex4::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{   
    ComplexType Value = Chi(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);

    if((ChiBit1 == ChiBit3 && MatsubaraNumber1 == MatsubaraNumber3) != 
       (ChiBit2 == ChiBit3 && MatsubaraNumber2 == MatsubaraNumber3)){
        ComplexType gg = g1(MatsubaraNumber1)*g2(MatsubaraNumber2);

        if(ChiBit1 == ChiBit3)
            Value -= Chi.getBeta() * gg;
        else
            Value += Chi.getBeta() * gg;
    }

    return Value;
}

ComplexType Vertex4::getAmputated(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    ComplexType Value = Chi(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
    ComplexType g1Value = g1(MatsubaraNumber1);
    ComplexType g2Value = g2(MatsubaraNumber2);
    ComplexType g3Value = g3(MatsubaraNumber3);
    ComplexType g4Value = g4(MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3);

    if((ChiBit1 == ChiBit3 && MatsubaraNumber1 == MatsubaraNumber3) != 
       (ChiBit2 == ChiBit3 && MatsubaraNumber2 == MatsubaraNumber3)){

        if(ChiBit1 == ChiBit3)
            Value -= Chi.getBeta() * g1Value*g2Value;
        else
            Value += Chi.getBeta() * g1Value*g2Value;
    }

    Value /= g1Value*g2Value*g3Value*g4Value;

    return Value;  
}