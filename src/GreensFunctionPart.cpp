#include "GreensFunctionPart.h"

GreensFunctionPart::greenTerm::greenTerm(ComplexType Residue, ComplexType Pole) : Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::greenTerm::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

GreensFunctionPart::GreensFunctionPart( AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                                        HamiltonianPart& HpartInner, HamiltonianPart& HpartOuter,
                                        DensityMatrixPart& DMpartInner, DensityMatrixPart& DMpartOuter) :
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        C(C), CX(CX)
{}

void GreensFunctionPart::compute(void)
{
    SparseMatrixType& Cmatrix = C.value();
    SparseMatrixType& CXmatrix = CX.value();
       
    QuantumState index1ket;
    QuantumState index1ketMax = CXmatrix.outerSize();
    
    for(index1ket=0; index1ket<index1ketMax; ++index1ket){
        SparseMatrixType::InnerIterator index2bra(CXmatrix,index1ket);       
        while(index2bra){
            QuantumState index2ket = index2bra.index();
            SparseMatrixType::InnerIterator index1bra(Cmatrix,index2ket);
            while(index1bra){
                //DEBUG("<" << index1bra.index() << "|c|" << 
                //      index2ket << "><" << index2bra.index() << "|c^+|" << index1ket << ">")
                QuantumState index1 = index1bra.index();
                if(index1 == index1ket){
                    ComplexType Residue = index1bra.value() * index2bra.value() * 
                                          (DMpartOuter.weight(index1) + DMpartInner.weight(index2ket));
                    ComplexType Pole = HpartInner.reV(index2ket) - HpartOuter.reV(index1);
              
                    Terms.push_back(greenTerm(Residue,Pole));
                }
                ++index1bra;
            }
            ++index2bra;
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
