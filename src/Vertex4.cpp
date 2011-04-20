#include "Vertex4.h"
#include <Eigen/LU> 

Vertex4::Vertex4(const BitClassification &IndexInfo, TwoParticleGFContainer &Chi, GFContainer &g) :
Chi(Chi), g(g), IndexInfo(IndexInfo), InvertedGFs(2*Chi.getNumberOfMatsubaras()+1), Storage(new TwoParticleGFPart::MatsubaraContainer(Chi.getBeta()))
{
    ParticleIndex N_bit = IndexInfo.getBitSize();
    NumberOfMatsubaras = Chi.getNumberOfMatsubaras();
    Storage->prepare(NumberOfMatsubaras);
    
    for(long MatsubaraNumber=-NumberOfMatsubaras; MatsubaraNumber<=NumberOfMatsubaras; ++MatsubaraNumber){
        MatrixType GMatrix(N_bit,N_bit);
        for(ParticleIndex index1=0; index1 < N_bit; ++index1)
        for(ParticleIndex index2=0; index2 < N_bit; ++index2){
            GMatrix(index1,index2) = g(index1,index2,MatsubaraNumber);
        }
        //DEBUG("InvertedGFs[" << MatsubaraNumber << "] = " << GMatrix.inverse() );
        InvertedGFs[MatsubaraNumber+NumberOfMatsubaras] = GMatrix.inverse();
    }

}


void Vertex4::compute(const TwoParticleGFContainer::IndexCombination& in)
{
    for (long MatsubaraNumber1=-NumberOfMatsubaras; MatsubaraNumber1 <= NumberOfMatsubaras; ++MatsubaraNumber1)
        for (long MatsubaraNumber2=-NumberOfMatsubaras; MatsubaraNumber1 <= NumberOfMatsubaras; ++MatsubaraNumber1)
            for (long MatsubaraNumber3=-NumberOfMatsubaras; MatsubaraNumber1 <= NumberOfMatsubaras; ++MatsubaraNumber1){
    
                ComplexType Value = Chi(in, MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
                RealType beta = Chi.getBeta();

                if(MatsubaraNumber1 == MatsubaraNumber3)
                    Value += beta*  g(in.Indices[0],in.Indices[2],MatsubaraNumber1)*
                                    g(in.Indices[1],in.Indices[3],MatsubaraNumber2);
                if(MatsubaraNumber1 == MatsubaraNumber2)
                    Value -= beta*  g(in.Indices[0],in.Indices[3],MatsubaraNumber1)*
                                    g(in.Indices[1],in.Indices[2],MatsubaraNumber2);

                (*Storage)(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3) = Value;
            }
}


ComplexType Vertex4::operator()(const TwoParticleGFContainer::IndexCombination& in,
                                long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    if(MatsubaraNumber1 > NumberOfMatsubaras || MatsubaraNumber1 < -NumberOfMatsubaras ||
       MatsubaraNumber2 > NumberOfMatsubaras || MatsubaraNumber2 < -NumberOfMatsubaras ||
       MatsubaraNumber3 > NumberOfMatsubaras || MatsubaraNumber3 < -NumberOfMatsubaras
       ) return 0;
    return (*Storage)(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3); 

    ComplexType Value = Chi(in, MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
    RealType beta = Chi.getBeta();

    if(MatsubaraNumber1 == MatsubaraNumber3)
        Value += beta*  g(in.Indices[0],in.Indices[2],MatsubaraNumber1)*
                        g(in.Indices[1],in.Indices[3],MatsubaraNumber2);
    if(MatsubaraNumber1 == MatsubaraNumber2)
        Value -= beta*  g(in.Indices[0],in.Indices[3],MatsubaraNumber1)*
                        g(in.Indices[1],in.Indices[2],MatsubaraNumber2);

    return Value;
}

ComplexType Vertex4::getAmputated(const TwoParticleGFContainer::IndexCombination& in,
                                  long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
   
    if(MatsubaraNumber1 > NumberOfMatsubaras || MatsubaraNumber1 < -NumberOfMatsubaras ||
       MatsubaraNumber2 > NumberOfMatsubaras || MatsubaraNumber2 < -NumberOfMatsubaras ||
       MatsubaraNumber3 > NumberOfMatsubaras || MatsubaraNumber3 < -NumberOfMatsubaras 
       ) return 0;

    ParticleIndex N_bit = IndexInfo.getBitSize();
    ComplexType Value = 0;
    
    for(ParticleIndex InnerIndex1=0; InnerIndex1 < N_bit; ++InnerIndex1)
    for(ParticleIndex InnerIndex2=0; InnerIndex2 < N_bit; ++InnerIndex2)
    for(ParticleIndex InnerIndex3=0; InnerIndex3 < N_bit; ++InnerIndex3)
    for(ParticleIndex InnerIndex4=0; InnerIndex4 < N_bit; ++InnerIndex4){
        TwoParticleGFContainer::IndexCombination iin(InnerIndex1,InnerIndex2,InnerIndex3,InnerIndex4);
        Value += (*this)(iin,MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3) *
                 InvertedGFs[MatsubaraNumber1+NumberOfMatsubaras](in.Indices[0],iin.Indices[0])*
                 InvertedGFs[MatsubaraNumber2+NumberOfMatsubaras](in.Indices[1],iin.Indices[1])*
                 InvertedGFs[MatsubaraNumber3+NumberOfMatsubaras](iin.Indices[2],in.Indices[2])*
                 InvertedGFs[MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3+NumberOfMatsubaras](iin.Indices[3],in.Indices[3]);
    }

    return Value;  
}
