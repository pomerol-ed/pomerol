#include "GFContainer.h"
extern ostream& OUTPUT_STREAM;

GFContainer::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cdagindex2)
{
    Indices[0]=cindex1;
    Indices[1]=cdagindex2;
}

bool GFContainer::IndexCombination::operator<(const GFContainer::IndexCombination& rhs) const
{
  return (Indices[0]<rhs.Indices[0]) || (Indices[0]==rhs.Indices[0] && Indices[1] < rhs.Indices[1]);
}


std::ostream& operator<<(std::ostream& output,const GFContainer::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}


/*=========================================================================*/

GFContainer::GFContainer ( StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,BitClassification& IndexInfo, FieldOperatorContainer& Operators):S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators)
{
};

MatrixType& GFContainer::operator()(long MatsubaraNumber)
{
    MatrixType* Output = new MatrixType(IndexInfo.getBitSize(), IndexInfo.getBitSize());
    Output->setZero();
    for (std::map<IndexCombination, GreensFunction*>::iterator it1=mapGreensFunctions.begin(); it1!=mapGreensFunctions.end(); it1++){
        (*Output)(it1->first.Indices[0],it1->first.Indices[1]) = (*it1->second)(MatsubaraNumber);
        };
    return *Output; 
};
