#include "Vertex4Part.h"

Vertex4Part::Vertex4TermType1::Vertex4TermType1(RealType weight, 
                                                RealType E1, RealType E2, RealType E3, RealType E4,
                                                Permutation3& Permutation)
{
    Residue = -weight*Permutation.sign;
    Poles[Permutation.perm[0]] = E1 - E4;
    Poles[Permutation.perm[1]] = E2 - E1;
    Poles[Permutation.perm[2]] = E3 - E2;
}

ComplexType Vertex4Part::Vertex4TermType1::operator()
            (ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const
{
    return Residue/((Frequency1 - Poles[0])*(Frequency2 - Poles[1])*(-Frequency3 - Poles[2]));
}

Vertex4Part::Vertex4Part(
                FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation) :
O1(O1), O2(O2), O3(O3), CX4(CX4), 
Hpart1(Hpart1), Hpart2(Hpart2), Hpart3(Hpart3), Hpart4(Hpart4),
DMpart1(DMpart1), DMpart2(DMpart2), DMpart3(DMpart3), DMpart4(DMpart4),
Permutation(Permutation)
{
    // TODO
}

void Vertex4Part::compute(void)
{
    SparseMatrixType& O1matrix = O1.value();
    SparseMatrixType& O2matrix = O2.value();    
    SparseMatrixType& O3matrix = O3.value();
    SparseMatrixType& CX4matrix = CX4.value();
    
    QuantumState index4ket;
    QuantumState index4ketMax = CX4matrix.outerSize();
    
    for(index4ket=0; index4ket<index4ketMax; ++index4ket){
        SparseMatrixType::InnerIterator index3bra(CX4matrix,index4ket);       
        while(index3bra){
            QuantumState index3ket = index3bra.index();
            SparseMatrixType::InnerIterator index2bra(O3matrix,index3ket);
            while(index2bra){
                QuantumState index2ket = index2bra.index();
                SparseMatrixType::InnerIterator index1bra(O2matrix,index2ket);
                while(index1bra){
                    QuantumState index1ket = index1bra.index();
                    SparseMatrixType::InnerIterator index4bra(O1matrix,index1ket);
                    while(index4bra){
                        QuantumState index4 = index4bra.index();
                        if(index4 == index4ket){
                            // TODO
                        }
                        ++index4bra;
                    }
                    ++index1bra;
                }
                ++index2bra;
            }
            ++index3bra;
        }
    }
}

ComplexType Vertex4Part::operator()(ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3)
{
    ComplexType Value = 0;
    for(std::list<Vertex4TermType1>::const_iterator pTerm = TermsType1.begin(); pTerm != TermsType1.end(); ++pTerm)
        Value += (*pTerm)(Frequency1,Frequency2,Frequency3);
    
    return Value;
}