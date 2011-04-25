#include "Vertex4.h"
#include <Eigen/LU> 

Vertex4::Vertex4(const IndexClassification &IndexInfo, TwoParticleGFContainer &Chi, GFContainer &g) :
Chi(Chi), g(g), IndexInfo(IndexInfo), InvertedGFs(2*Chi.getNumberOfMatsubaras()+1)
{
    NumberOfMatsubaras = Chi.getNumberOfMatsubaras();
}

ComplexType Vertex4::operator()(const IndexCombination& in,
                    long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    return this->getAmputatedValue(in,MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}
//============================= UnAmputated methods ==============================//

void Vertex4::prepareUnAmputated()
{
    INFO("Preparing unamputated value container map ...");

    for (std::vector<IndexCombination*>::const_iterator it1=Chi.getNonTrivialCombinations().begin(); it1!=Chi.getNonTrivialCombinations().end(); ++it1){
        mapUnAmputatedValues[**it1]= new TwoParticleGFPart::MatsubaraContainer(Chi.getBeta());
        mapUnAmputatedValues[**it1]->prepare(NumberOfMatsubaras);
    };
    for (std::vector<IndexCombination*>::const_iterator it1=Chi.getTrivialCombinations().begin(); it1!=Chi.getTrivialCombinations().end(); ++it1){
        mapUnAmputatedValues[**it1]= new TwoParticleGFPart::MatsubaraContainer(Chi.getBeta());
        mapUnAmputatedValues[**it1]->prepare(NumberOfMatsubaras);
    };
};

void Vertex4::computeUnAmputated(const IndexCombination& in)
{
    INFO("Computing unamputated values for " << in << " combination");
    for (long MatsubaraNumber1=-NumberOfMatsubaras; MatsubaraNumber1 < NumberOfMatsubaras; ++MatsubaraNumber1)
        for (long MatsubaraNumber2=-NumberOfMatsubaras; MatsubaraNumber2 < NumberOfMatsubaras; ++MatsubaraNumber2)
            for (long MatsubaraNumber3=-NumberOfMatsubaras; MatsubaraNumber3 < NumberOfMatsubaras; ++MatsubaraNumber3)
                if ( MatsubaraNumber1 + MatsubaraNumber2 - MatsubaraNumber3 < NumberOfMatsubaras && MatsubaraNumber1 + MatsubaraNumber2 - MatsubaraNumber3 >= -NumberOfMatsubaras){ 
                    ComplexType Value = Chi(in, MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
                    RealType beta = Chi.getBeta();

                    if(MatsubaraNumber1 == MatsubaraNumber3)
                        Value += beta*  g(in.Indices[0],in.Indices[2],MatsubaraNumber1)*
                                        g(in.Indices[1],in.Indices[3],MatsubaraNumber2);
                    if(MatsubaraNumber2 == MatsubaraNumber3)
                        Value -= beta*  g(in.Indices[0],in.Indices[3],MatsubaraNumber1)*
                                        g(in.Indices[1],in.Indices[2],MatsubaraNumber2);

                    (*mapUnAmputatedValues[in])(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3) = Value;
            }
}

void Vertex4::computeUnAmputated()
{
    for (std::map<IndexCombination,TwoParticleGFPart::MatsubaraContainer*>::iterator it=mapUnAmputatedValues.begin();
            it!=mapUnAmputatedValues.end(); ++it){
                computeUnAmputated(it->first);
                }
};


ComplexType Vertex4::getUnAmputatedValue(const IndexCombination& in,
                                long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    if(MatsubaraNumber1 > NumberOfMatsubaras || MatsubaraNumber1 < -NumberOfMatsubaras ||
       MatsubaraNumber2 > NumberOfMatsubaras || MatsubaraNumber2 < -NumberOfMatsubaras ||
       MatsubaraNumber3 > NumberOfMatsubaras || MatsubaraNumber3 < -NumberOfMatsubaras ||
       mapUnAmputatedValues.count(in) == 0
       ) return 0;
    return (*mapUnAmputatedValues[in])(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);
}

//============================= Amputated methods ==============================//

void Vertex4::prepareAmputated(std::vector<IndexCombination*>& in)
{
    NonTrivialAmputatedCombinations=in;
    INFO("Inverting Greens Function...");
    ParticleIndex N_bit = IndexInfo.getIndexSize();
    for(long MatsubaraNumber=-NumberOfMatsubaras; MatsubaraNumber<=NumberOfMatsubaras; ++MatsubaraNumber){
        MatrixType GMatrix(N_bit,N_bit);
        for(ParticleIndex index1=0; index1 < N_bit; ++index1)
        for(ParticleIndex index2=0; index2 < N_bit; ++index2){
            GMatrix(index1,index2) = g(index1,index2,MatsubaraNumber);
        }
        //DEBUG("InvertedGFs[" << MatsubaraNumber << "] = " << GMatrix.inverse() );
        InvertedGFs[MatsubaraNumber+NumberOfMatsubaras] = GMatrix.inverse();
    }
    INFO("Preparing amputated value container map ...");
    for (std::vector<IndexCombination*>::const_iterator it1=NonTrivialAmputatedCombinations.begin(); 
            it1!=NonTrivialAmputatedCombinations.end(); ++it1){
                mapAmputatedValues[**it1]= new TwoParticleGFPart::MatsubaraContainer(Chi.getBeta());
                mapAmputatedValues[**it1]->prepare(NumberOfMatsubaras);
    }       ;
};

void Vertex4::computeAmputated()
{
    for (std::map<IndexCombination,TwoParticleGFPart::MatsubaraContainer*>::iterator it=mapAmputatedValues.begin();
            it!=mapAmputatedValues.end(); ++it){
                    computeAmputated(it->first);
                    };
}

void Vertex4::computeAmputated(const IndexCombination& in)
{
    INFO("Computing amputated values for " << in << " combination");
    for (long MatsubaraNumber1=-NumberOfMatsubaras; MatsubaraNumber1 < NumberOfMatsubaras; ++MatsubaraNumber1)
        for (long MatsubaraNumber2=-NumberOfMatsubaras; MatsubaraNumber2 < NumberOfMatsubaras; ++MatsubaraNumber2)
            for (long MatsubaraNumber3=-NumberOfMatsubaras; MatsubaraNumber3 < NumberOfMatsubaras; ++MatsubaraNumber3)
                if ( MatsubaraNumber1 + MatsubaraNumber2 - MatsubaraNumber3 < NumberOfMatsubaras && MatsubaraNumber1 + MatsubaraNumber2 - MatsubaraNumber3 >= -NumberOfMatsubaras){ 
                    ComplexType Value=0;
                    for (std::map<IndexCombination,TwoParticleGFPart::MatsubaraContainer*>::iterator it=mapUnAmputatedValues.begin();
                        it!=mapUnAmputatedValues.end(); ++it){
                            IndexCombination iin = it->first;

                            Value += this->getUnAmputatedValue(iin,MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3) *
                                    InvertedGFs[MatsubaraNumber1+NumberOfMatsubaras](in.Indices[0],iin.Indices[0])*
                                    InvertedGFs[MatsubaraNumber2+NumberOfMatsubaras](in.Indices[1],iin.Indices[1])*
                                    InvertedGFs[MatsubaraNumber3+NumberOfMatsubaras](iin.Indices[2],in.Indices[2])*
                                    InvertedGFs[MatsubaraNumber1+MatsubaraNumber2-MatsubaraNumber3+NumberOfMatsubaras](iin.Indices[3],in.Indices[3]);
                            };
                    (*mapAmputatedValues[in])(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3) = Value;
                    }
}


ComplexType Vertex4::getAmputatedValue(const IndexCombination& in,
                                  long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
   
    if(MatsubaraNumber1 > NumberOfMatsubaras || MatsubaraNumber1 < -NumberOfMatsubaras ||
       MatsubaraNumber2 > NumberOfMatsubaras || MatsubaraNumber2 < -NumberOfMatsubaras ||
       MatsubaraNumber3 > NumberOfMatsubaras || MatsubaraNumber3 < -NumberOfMatsubaras || 
       mapAmputatedValues.count(in) == 0
       ) return 0;

    return (*mapAmputatedValues[in])(MatsubaraNumber1, MatsubaraNumber2, MatsubaraNumber3);
};


