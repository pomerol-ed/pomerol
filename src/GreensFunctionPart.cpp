#include "GreensFunctionPart.h"

GreensFunctionPart::GreensTerm::GreensTerm(ComplexType Residue, ComplexType Pole) : Residue(Residue), Pole(Pole) {};
ComplexType GreensFunctionPart::GreensTerm::operator()(ComplexType Frequency) const { return Residue/(Frequency - Pole); }

std::ostream& operator<<(std::ostream& out, const GreensFunctionPart::GreensTerm& Term)
{
    out << Term.Residue << "/(z - " << Term.Pole << ")";
    return out;
}

GreensFunctionPart::GreensFunctionPart( AnnihilationOperatorPart& C, CreationOperatorPart& CX, 
                                        HamiltonianPart& HpartInner, HamiltonianPart& HpartOuter,
                                        DensityMatrixPart& DMpartInner, DensityMatrixPart& DMpartOuter) :
                                        HpartInner(HpartInner), HpartOuter(HpartOuter),
                                        DMpartInner(DMpartInner), DMpartOuter(DMpartOuter),
                                        C(C), CX(CX)
{}
/*
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
                QuantumState index1 = index1bra.index();
                if(index1 == index1ket){
                    ComplexType Residue = index1bra.value() * index2bra.value() * 
                                          (DMpartOuter.weight(index1) + DMpartInner.weight(index2ket));
                    ComplexType Pole = HpartInner.reV(index2ket) - HpartOuter.reV(index1);
                    
                    if(abs(Residue) > MATRIX_ELEMENT_TOLERANCE)
                        Terms.push_back(GreensTerm(Residue,Pole));
                }
                ++index1bra;
            }
            ++index2bra;
        }
    }
}
*/
void GreensFunctionPart::compute(void)
{
    Terms.clear();
  
    RowMajorMatrixType& Cmatrix = C.getRowMajorValue();
    ColMajorMatrixType& CXmatrix = CX.getColMajorValue();
    QuantumState outerSize = Cmatrix.outerSize();
      
    for(QuantumState index1=0; index1<outerSize; ++index1){
        RowMajorMatrixType::InnerIterator Cinner(Cmatrix,index1);
        ColMajorMatrixType::InnerIterator CXinner(CXmatrix,index1);
        
        while(Cinner && CXinner){
            QuantumState C_index2 = Cinner.index();
            QuantumState CX_index2 = CXinner.index();

            if(C_index2 == CX_index2){
                ComplexType Residue = Cinner.value() * CXinner.value() * 
                                      (DMpartOuter.weight(index1) + DMpartInner.weight(C_index2));
                if(abs(Residue) > MATRIX_ELEMENT_TOLERANCE)
		{
                	ComplexType Pole = HpartInner.reV(C_index2) - HpartOuter.reV(index1);
                	Terms.push_back(GreensTerm(Residue,Pole));
		};
                ++Cinner;
                ++CXinner;
            }else{
                if(CX_index2 < C_index2) for(;QuantumState(CXinner.index())<C_index2; ++CXinner);
                else for(;QuantumState(Cinner.index())<CX_index2; ++Cinner);
            }
        }
    }
}



ComplexType GreensFunctionPart::operator()(long MatsubaraNumber) const
{
    ComplexType MatsubaraSpacing = I*M_PI/DMpartInner.getBeta();
  
    ComplexType G = 0;
    for(std::list<GreensTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
        G += (*pTerm)(MatsubaraSpacing*RealType(2*MatsubaraNumber+1));
        
    return G;
}

const std::list<GreensFunctionPart::GreensTerm>& GreensFunctionPart::getTerms(void) const
{
    return Terms;
}
