#include "GreensFunctionPart.h"

GreensFunctionPart::greenTerm::greenTerm(ComplexType Residue, ComplexType Pole) : Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::greenTerm::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

GreensFunctionPart::GreensFunctionPart( AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                                        HamiltonianPart& Hpart, DensityMatrixPart& DMpart) :
                                        Hpart(Hpart), C(C), CX(CX), DMpart(DMpart)
{}

void GreensFunctionPart::compute(void)
{
    RowMajorMatrixType& Cmatrix = C.value();
    ColMajorMatrixType& CXmatrix = CX.value();
    QuantumState outerSize = Cmatrix.outerSize();
      
    for(QuantumState m=0; m<outerSize; ++m){
        RowMajorMatrixType::InnerIterator Cinner(Cmatrix,m);
        ColMajorMatrixType::InnerIterator CXinner(CXmatrix,m);
        
        while(Cinner && CXinner){
            QuantumState C_n = Cinner.index();
            QuantumState CX_n = CXinner.index();
            if(C_n == CX_n){
                ComplexType Residue = Cinner.value() * CXinner.value() * 
                                      (DMpart.weight(m) + DMpart.weight(C_n));
                ComplexType Pole = Hpart.reV(C_n) - Hpart.reV(m);
              
                Terms.push_back(greenTerm(Residue,Pole));
              
                ++Cinner;
                ++CXinner;
            }else{
                if(CX_n < C_n) for(;QuantumState(CXinner.index())<C_n; ++CXinner);
                else for(;QuantumState(Cinner.index())<CX_n; ++Cinner);
            }
        }
    }
}

ComplexType GreensFunctionPart::operator()(ComplexType Frequency)
{
    ComplexType G = 0;
    for(std::list<greenTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
        G += (*pTerm)(Frequency);
        
    return G;
}
