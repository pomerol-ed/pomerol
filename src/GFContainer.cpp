#include "GFContainer.h"
extern std::ostream& OUTPUT_STREAM;

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
output << "(" << out.Indices[0] << out.Indices[1] << ")" << std::flush;
return output;
}


/*=========================================================================*/

GFContainer::GFContainer ( StatesClassification &S, Hamiltonian &H, DensityMatrix &DM,IndexClassification& IndexInfo, FieldOperatorContainer& Operators):ComputableObject(),S(S),H(H),DM(DM),IndexInfo(IndexInfo),Operators(Operators)
{
};

MatrixType& GFContainer::operator()(long MatsubaraNumber)
{
    MatrixType* Output = new MatrixType(IndexInfo.getIndexSize(), IndexInfo.getIndexSize());
    Output->setZero();
    for (std::map<IndexCombination, GreensFunction*>::iterator it1=mapGreensFunctions.begin(); it1!=mapGreensFunctions.end(); it1++){
        (*Output)(it1->first.Indices[0],it1->first.Indices[1]) = (*it1->second)(MatsubaraNumber);
        };
    return *Output; 
};

ComplexType GFContainer::operator()(ParticleIndex i, ParticleIndex j, long MatsubaraNumber)
{
    return (!vanishes(i,j))?(*mapGreensFunctions[IndexCombination(i,j)])(MatsubaraNumber):0;
};

bool GFContainer::vanishes(ParticleIndex i, ParticleIndex j)
{
   return (mapGreensFunctions.count(IndexCombination(i,j))==0); 
};

void GFContainer::prepare()
{
if (Status == Constructed){
    for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); ++i)
        for (ParticleIndex j=0; j<IndexInfo.getIndexSize(); ++j){
                GreensFunction *GF = new GreensFunction (S,H,Operators.getAnnihilationOperator(i),Operators.getCreationOperator(j),DM);
                GF->prepare();
                if (!GF->vanishes()) mapGreensFunctions[IndexCombination(i,j)]=GF;
            }
    Status = Prepared;
    };
};

void GFContainer::compute()
{
if (Status == Prepared){
    for (std::map<IndexCombination,GreensFunction*>::iterator it1=mapGreensFunctions.begin();it1!=mapGreensFunctions.end();++it1){
           DEBUG("GFContainer: computing G_{" << it1->first << "}");
           it1->second->compute(); 
        };
    Status = Computed;
    };
};

void GFContainer::dumpToPlainText(long wn)
{
    for (std::map<IndexCombination,GreensFunction*>::iterator it1=mapGreensFunctions.begin();it1!=mapGreensFunctions.end();++it1){
           it1->second->dumpToPlainText(wn); 
        };
};
