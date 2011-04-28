#include "FourIndexObject.h"
//
//FourIndexObject::IndexCombination
//

FourIndexObject::IndexCombination::IndexCombination(ParticleIndex cindex1, ParticleIndex cindex2, ParticleIndex cdagindex3, ParticleIndex cdagindex4)
{
    Indices[0]=cindex1;
    Indices[1]=cindex2;
    Indices[2]=cdagindex3;
    Indices[3]=cdagindex4;
}

bool FourIndexObject::IndexCombination::operator<(const FourIndexObject::IndexCombination& rhs) const
{
  return (Indices[0] < rhs.Indices[0]) || 
         (Indices[0] == rhs.Indices[0] && Indices[1] < rhs.Indices[1] ) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] < rhs.Indices[2]) ||
         (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] < rhs.Indices[3]); 
}

bool FourIndexObject::IndexCombination::operator==(const FourIndexObject::IndexCombination& rhs) const
{
    return (Indices[0] == rhs.Indices[0] && Indices[1] == rhs.Indices[1] && Indices[2] == rhs.Indices[2] && Indices[3] == rhs.Indices[3]); 
}

bool FourIndexObject::IndexCombination::operator!=(const FourIndexObject::IndexCombination& rhs) const
{
    return !(*this==rhs);
}

std::ostream& operator<<(std::ostream& output,const FourIndexObject::IndexCombination& out)
{
output << "(" << out.Indices[0] << out.Indices[1] << out.Indices[2] << out.Indices[3] << ")";
return output;
}

//
// Matsubara Container
//
FourIndexObject::MatsubaraContainer::MatsubaraContainer(RealType beta):MatsubaraSpacing(I*M_PI/beta){};

void FourIndexObject::MatsubaraContainer::prepare(long NumberOfMatsubaras_)
{
    NumberOfMatsubaras=NumberOfMatsubaras_;
    Data.resize(4*NumberOfMatsubaras);
    FermionicFirstIndex.resize(4*NumberOfMatsubaras);
    for (int BosonicIndex=-2*NumberOfMatsubaras;BosonicIndex<=(int)(2*NumberOfMatsubaras)-2;BosonicIndex++)
    {
        int Size=(BosonicIndex+1>0)?BosonicIndex+1-NumberOfMatsubaras:-NumberOfMatsubaras;
        FermionicFirstIndex[BosonicIndex+2*NumberOfMatsubaras]=Size;
        Size=((BosonicIndex+1<0)?BosonicIndex+NumberOfMatsubaras:NumberOfMatsubaras-1) - Size + 1;
        Size=(Size<=0)?0:Size;
        //      cout << "Freq = " << BosonicIndex << ", Size = " << Size << ", First Index = " << FermionicFirstIndex[BosonicIndex+2*NumberOfMatsubaras] << endl;
        Data[BosonicIndex+2*NumberOfMatsubaras].resize(Size,Size);
        Data[BosonicIndex+2*NumberOfMatsubaras].setZero();
    };
    //    exit(0);
};

FourIndexObject::MatsubaraContainer& FourIndexObject::MatsubaraContainer::operator+= (const MatsubaraContainer& rhs)
{
    for (long BosonicIndex=0;BosonicIndex<=(4*NumberOfMatsubaras)-2;BosonicIndex++){
        Data[BosonicIndex]+=rhs.Data[BosonicIndex];
    }
    return (*this);
};

void FourIndexObject::MatsubaraContainer::clear()
{
    for (long BosonicIndex=0;BosonicIndex<=(4*NumberOfMatsubaras)-2;BosonicIndex++){
        Data[BosonicIndex].resize(0,0);
    }
}

void FourIndexObject::MatsubaraContainer::fill(const std::list<TwoParticleGFPart::NonResonantTerm>& NonResonantTerms, const std::list<TwoParticleGFPart::ResonantTerm>& ResonantTerms, Permutation3 Permutation)
{
    for (long BosonicIndex=0;BosonicIndex<=(4*NumberOfMatsubaras)-2;BosonicIndex++){
        for (long nuIndex=0;nuIndex<Data[BosonicIndex].cols();++nuIndex){
            for (long nu1Index=0;nu1Index<Data[BosonicIndex].cols();++nu1Index){
                
                long FermionicIndexShift = FermionicFirstIndex[BosonicIndex];
                long MatsubaraNumber2 = nuIndex +FermionicIndexShift;
                long MatsubaraNumber1 = BosonicIndex-2*NumberOfMatsubaras-MatsubaraNumber2;
                long MatsubaraNumber3 = nu1Index+FermionicIndexShift;
                
                long MatsubaraNumberOdd1 = 2*MatsubaraNumber1 + 1;
                long MatsubaraNumberOdd2 = 2*MatsubaraNumber2 + 1;
                long MatsubaraNumberOdd3 = 2*MatsubaraNumber3 + 1;
                ComplexType Frequencies[3] = {  MatsubaraSpacing * RealType(MatsubaraNumberOdd1),
                                                MatsubaraSpacing * RealType(MatsubaraNumberOdd2),
                                               -MatsubaraSpacing * RealType(MatsubaraNumberOdd3)};
                                    
                ComplexType z1 = Frequencies[Permutation.perm[0]];                                    
                ComplexType z2 = Frequencies[Permutation.perm[1]];
                ComplexType z3 = Frequencies[Permutation.perm[2]];
    
                ComplexType Value = 0;
                for(std::list<TwoParticleGFPart::NonResonantTerm>::const_iterator pTerm = NonResonantTerms.begin(); pTerm != NonResonantTerms.end(); ++pTerm)
                    Value += (*pTerm)(z1,z2,z3);
                for(std::list<TwoParticleGFPart::ResonantTerm>::const_iterator pTerm = ResonantTerms.begin(); pTerm != ResonantTerms.end(); ++pTerm)
                    Value += (*pTerm)(z1,z2,z3);

                Data[BosonicIndex](nuIndex,nu1Index) += Value;
            };
        };
    };
};


//
// FourIndexSingleObject
//

//
// FourIndexContainerObject
//
const Permutation4 FourIndexContainerObject::TrivialOperatorPermutations[4] = { permutations4[0], permutations4[1], permutations4[6], permutations4[7] };
