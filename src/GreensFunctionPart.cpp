#include "GreensFunctionPart.h"

GreensFunctionPart::greenTerm::greenTerm(ComplexType Residue, ComplexType Pole) : Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::greenTerm::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

GreensFunctionPart::GreensFunctionPart( AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                                        HamiltonianPart& Hpart, DensityMatrixPart& DMpart)
{
    // TODO
}

ComplexType GreensFunctionPart::operator()(ComplexType Frequency)
{
    ComplexType G = 0;
    for(std::list<greenTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
        G += (*pTerm)(Frequency);
        
    return G;
}
