#include "TwoParticleGFPart.h"

//DEBUG
long term_counter=0;

// Make the lagging index catch up or outrun the leading index.
inline bool chaseIndices(RowMajorMatrixType::InnerIterator& index1_iter, 
                         ColMajorMatrixType::InnerIterator& index2_iter)
{
    InnerQuantumState index1 = index1_iter.index();
    InnerQuantumState index2 = index2_iter.index();

    if(index1 == index2) return true;

    if(index1 < index2) 
        for(;InnerQuantumState(index1_iter.index())<index2 && index1_iter; ++index1_iter);
    else 
        for(;InnerQuantumState(index2_iter.index())<index1 && index2_iter; ++index2_iter);

    return false;
}

//
// TwoParticleGFPart::TwoParticleGFTerm
//
inline
TwoParticleGFPart::TwoParticleGFTerm::TwoParticleGFTerm(ComplexType Coeff, RealType beta,
                                RealType Ei, RealType Ej, RealType Ek, RealType El,
                                RealType Wi, RealType Wj, RealType Wk, RealType Wl,
                                Permutation3& Permutation)
{
    z1 = Permutation.perm[0];
    z2 = Permutation.perm[1];
    z3 = Permutation.perm[2];

    Poles[0] = Ej - Ei;
    Poles[1] = Ek - Ej;
    Poles[2] = El - Ek;

    CoeffZ2 = -Coeff*(Wj + Wk);
    CoeffZ4 = Coeff*(Wi + Wl);
    CoeffZ1Z2Res = Coeff*beta*Wi;
    CoeffZ1Z2NonRes = Coeff*(Wk - Wi);
    CoeffZ2Z3Res = -Coeff*beta*Wj;
    CoeffZ2Z3NonRes = Coeff*(Wj - Wl);
}

ComplexType TwoParticleGFPart::TwoParticleGFTerm::operator()(
                                ComplexType Frequency1, ComplexType Frequency2, ComplexType Frequency3) const
{
    ComplexType w[3] = {Frequency1,Frequency2,-Frequency3};

    ComplexType Value = 0;
    ComplexType Z1MinusP1 = w[z1] - Poles[0];
    ComplexType Z2MinusP2 = w[z2] - Poles[1];
    ComplexType Z3MinusP3 = w[z3] - Poles[2];

    Value += CoeffZ2/Z2MinusP2;
    Value += CoeffZ4/(Z1MinusP1+Z2MinusP2+Z3MinusP3);

    ComplexType Z1Z2P1P2 = Z1MinusP1 + Z2MinusP2;
    ComplexType Z2Z3P2P3 = Z2MinusP2 + Z3MinusP3;
    // 2 cases for each of '12' and '23' pairs: resonant and non-resonant.
    Value += (abs(Z1Z2P1P2) < ResonanceTolerance) ? CoeffZ1Z2Res : CoeffZ1Z2NonRes/Z1Z2P1P2;
    Value += (abs(Z2Z3P2P3) < ResonanceTolerance) ? CoeffZ2Z3Res : CoeffZ2Z3NonRes/Z2Z3P2P3;

    Value /= (Z1MinusP1*Z3MinusP3);
    return Value;
}

inline
bool TwoParticleGFPart::TwoParticleGFTerm::IsRelevant(const ComplexType &MatrixElementProd)
{
    return abs(MatrixElementProd) > MatrixElementTolerance;
}

//
// Matsubara Container
//
TwoParticleGFPart::MatsubaraContainer::MatsubaraContainer(RealType beta):MatsubaraSpacing(I*M_PI/beta){};

void TwoParticleGFPart::MatsubaraContainer::prepare(long NumberOfMatsubaras_)
{
    NumberOfMatsubaras=NumberOfMatsubaras_;
    Data.resize(4*NumberOfMatsubaras);
    FermionicFirstIndex.resize(4*NumberOfMatsubaras);
    for (int BosonicIndex=-2*NumberOfMatsubaras;BosonicIndex<=(int)(2*NumberOfMatsubaras)-2;BosonicIndex++)
    {
      int Size=(BosonicIndex+1-NumberOfMatsubaras>-NumberOfMatsubaras)?BosonicIndex+1-NumberOfMatsubaras:-NumberOfMatsubaras;
      FermionicFirstIndex[BosonicIndex+2*NumberOfMatsubaras]=Size;
      Size=((BosonicIndex+NumberOfMatsubaras<NumberOfMatsubaras-1)?BosonicIndex+NumberOfMatsubaras:NumberOfMatsubaras-1) - Size + 1;
      Size=(Size<=0)?0:Size;
//      cout << "Freq = " << BosonicIndex << ", Size = " << Size << ", First Index = " << FermionicFirstIndex[BosonicIndex+2*NumberOfMatsubaras] << endl;
      Data[BosonicIndex+2*NumberOfMatsubaras].resize(Size,Size);
      Data[BosonicIndex+2*NumberOfMatsubaras].setZero();
    };
//    exit(0);
};

ComplexType& TwoParticleGFPart::MatsubaraContainer::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
// {OMEGA,nu,nu' : OMEGA=w1+w2, nu=w1, nu'=w4=OMEGA-w3
{
    unsigned int RealBosonicIndex = MatsubaraNumber1 + MatsubaraNumber2 + 2*NumberOfMatsubaras;
    int nuIndex = MatsubaraNumber1-FermionicFirstIndex[RealBosonicIndex];
    int nu1Index= RealBosonicIndex-2*NumberOfMatsubaras-MatsubaraNumber3-FermionicFirstIndex[RealBosonicIndex];
    //cout << "Bosonic index : " << RealBosonicIndex - 2*NumberOfMatsubaras<< " shift : " << FermionicFirstIndex[RealBosonicIndex] << endl;
    if (nuIndex >= 0 && nuIndex < Data[RealBosonicIndex].rows() && nu1Index >= 0 && nu1Index < Data[RealBosonicIndex].rows() )
        return Data[RealBosonicIndex](nuIndex,nu1Index);
    else {
        static ComplexType temp(0.0,0.0);
        cout << "Warning! Index (" << MatsubaraNumber1 << "," << MatsubaraNumber2 << "," << MatsubaraNumber3 << "," << MatsubaraNumber1+ MatsubaraNumber2 - MatsubaraNumber3 << ") of TwoParticleGFPart is out of range, returning 0" << endl;
        return temp;
         };
};

TwoParticleGFPart::MatsubaraContainer& TwoParticleGFPart::MatsubaraContainer::operator+= (const MatsubaraContainer& rhs)
{
    for (long BosonicIndex=0;BosonicIndex<=(4*NumberOfMatsubaras)-2;BosonicIndex++){
        Data[BosonicIndex]+=rhs.Data[BosonicIndex];
    }
    return (*this);
};

inline void TwoParticleGFPart::MatsubaraContainer::fill(std::list<TwoParticleGFTerm> &Terms)
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
                ComplexType Frequency1 = MatsubaraSpacing * RealType(MatsubaraNumberOdd1);
                ComplexType Frequency2 = MatsubaraSpacing * RealType(MatsubaraNumberOdd2);
                ComplexType Frequency3 = MatsubaraSpacing * RealType(MatsubaraNumberOdd3);

                ComplexType Value = 0;
                for(std::list<TwoParticleGFTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
                    Value += (*pTerm)(Frequency1,Frequency2,Frequency3);

                Data[BosonicIndex](nuIndex,nu1Index)+=Value;
            };
        };
    };
};


void TwoParticleGFPart::MatsubaraContainer::clear()
{
    for (long BosonicIndex=0;BosonicIndex<=(4*NumberOfMatsubaras)-2;BosonicIndex++){
        Data[BosonicIndex].resize(0,0);
    }
}
//
// TwoParticleGFPart
//
TwoParticleGFPart::TwoParticleGFPart(
                FieldOperatorPart& O1, FieldOperatorPart& O2, FieldOperatorPart& O3, CreationOperatorPart& CX4,
                HamiltonianPart& Hpart1, HamiltonianPart& Hpart2, HamiltonianPart& Hpart3, HamiltonianPart& Hpart4,
                DensityMatrixPart& DMpart1, DensityMatrixPart& DMpart2, DensityMatrixPart& DMpart3, DensityMatrixPart& DMpart4,
                Permutation3 Permutation) :
O1(O1), O2(O2), O3(O3), CX4(CX4), 
Hpart1(Hpart1), Hpart2(Hpart2), Hpart3(Hpart3), Hpart4(Hpart4),
DMpart1(DMpart1), DMpart2(DMpart2), DMpart3(DMpart3), DMpart4(DMpart4),
Permutation(Permutation),Storage(new MatsubaraContainer(DMpart1.getBeta())){};

void TwoParticleGFPart::compute(long NumberOfMatsubaras, TwoParticleGFPart::ComputationMethod method)
{
    Terms.clear();

    switch(method){
        case ChasingIndices1: computeChasing1(NumberOfMatsubaras); break;
        case ChasingIndices2: computeChasing2(NumberOfMatsubaras); break;
        default: assert(0);
    };
};

void TwoParticleGFPart::clear()
{
    Terms.clear();
    Storage->clear();
};

void TwoParticleGFPart::computeChasing1(long NumberOfMatsubaras)
    // I don't have any pen now, so I'm writing here:
    // <1 | O1 | 2> <2 | O2 | 3> <3 | O3 |4> <4| CX4 |1>
    // Iterate over all values of |1><1|
    // Chase indices |3> and <3|
{
    Storage->prepare(NumberOfMatsubaras);
    RealType beta = DMpart1.getBeta();  

    RowMajorMatrixType& O1matrix = O1.getRowMajorValue();
    RowMajorMatrixType& O2matrix = O2.getRowMajorValue();    
    ColMajorMatrixType& O3matrix = O3.getColMajorValue();
    ColMajorMatrixType& CX4matrix = CX4.getColMajorValue();

    InnerQuantumState index1;
    InnerQuantumState index1Max = CX4matrix.outerSize();

    for(index1=0; index1<index1Max; ++index1){
        ColMajorMatrixType::InnerIterator index4bra(CX4matrix,index1);

        RealType E1 = Hpart1.reV(index1);
        RealType weight1 = DMpart1.weight(index1); 

        while(index4bra){
            InnerQuantumState index4ket = index4bra.index();
            RowMajorMatrixType::InnerIterator index2ket(O1matrix,index1);

            RealType E4 = Hpart4.reV(index4ket);
            RealType weight4 = DMpart4.weight(index4ket);

            while (index2ket){
                InnerQuantumState index2bra = index2ket.index();

                RealType E2 = Hpart2.reV(index2bra);
                RealType weight2 = DMpart2.weight(index2bra);

                RowMajorMatrixType::InnerIterator index3ket_iter(O2matrix,index2bra);
                ColMajorMatrixType::InnerIterator index3bra_iter(O3matrix,index4ket);
                while(index3bra_iter && index3ket_iter){
                    if(chaseIndices(index3ket_iter,index3bra_iter)){

                        InnerQuantumState index3 = index3ket_iter.index();

                        RealType E3 = Hpart3.reV(index3);                       
                        RealType weight3 = DMpart3.weight(index3);

                        ComplexType MatrixElement = index2ket.value()*
                                                    index3ket_iter.value()*
                                                    index3bra_iter.value()*
                                                    index4bra.value();

                        MatrixElement *= Permutation.sign;

                        if(TwoParticleGFTerm::IsRelevant(MatrixElement)){
                            term_counter++;
                            Terms.push_back(TwoParticleGFTerm(MatrixElement,beta,
                                                              E1,E2,E3,E4,
                                                              weight1,weight2,weight3,weight4,Permutation));
                        }

                        ++index3ket_iter;
                        ++index3bra_iter;
                    }
                };
                 ++index2ket;
            };
            ++index4bra;
        };
    };
    Storage->fill(Terms);
};

void TwoParticleGFPart::computeChasing2(long NumberOfMatsubaras)
{
    RealType beta = DMpart1.getBeta();
    Storage->prepare(NumberOfMatsubaras);
    // I don't have any pen now, so I'm writing here:
    // <1 | O1 | 2> <2 | O2 | 3> <3 | O3 |4> <4| CX4 |1>
    // Iterate over all values of |1><1| and |3><3|
    // Chase indices |2> and <2|, |4> and <4|.
    RowMajorMatrixType& O1matrix = O1.getRowMajorValue();
    ColMajorMatrixType& O2matrix = O2.getColMajorValue();    
    RowMajorMatrixType& O3matrix = O3.getRowMajorValue();
    ColMajorMatrixType& CX4matrix = CX4.getColMajorValue();

    InnerQuantumState index1;
//    #ifndef DMTruncate
    InnerQuantumState index1Max = CX4matrix.outerSize();
//    #else
//    InnerQuantumState index1Max = DMpart1.getMaximumTruncationState();
//    #endif


    InnerQuantumState index3;
    InnerQuantumState index3Max = O2matrix.outerSize();

    std::list<InnerQuantumState> Index4List;

    for(index1=0; index1<index1Max; ++index1){
    for(index3=0; index3<index3Max; ++index3){
        ColMajorMatrixType::InnerIterator index4bra_iter(CX4matrix,index1);       
        RowMajorMatrixType::InnerIterator index4ket_iter(O3matrix,index3);
        Index4List.clear();
        while (index4bra_iter && index4ket_iter){
            if(chaseIndices(index4ket_iter,index4bra_iter)){
                Index4List.push_back(index4bra_iter.index());
                ++index4bra_iter;
                ++index4ket_iter;
            }
        };

        if (!Index4List.empty())
        {
            RealType E1 = Hpart1.reV(index1);
            RealType E3 = Hpart3.reV(index3);
            RealType weight1 = DMpart1.weight(index1);
            RealType weight3 = DMpart3.weight(index3);

            ColMajorMatrixType::InnerIterator index2bra_iter(O2matrix,index3);
            RowMajorMatrixType::InnerIterator index2ket_iter(O1matrix,index1);       
            while (index2bra_iter && index2ket_iter){
                if (chaseIndices(index2ket_iter,index2bra_iter)){

                    InnerQuantumState index2 = index2ket_iter.index();
                    RealType E2 = Hpart2.reV(index2);
                    RealType weight2 = DMpart2.weight(index2);

                    for (std::list<InnerQuantumState>::iterator pIndex4 = Index4List.begin(); pIndex4!=Index4List.end(); ++pIndex4) 
                    {
                        InnerQuantumState index4 = *pIndex4;
                        RealType E4 = Hpart4.reV(index4);                       
                        RealType weight4 = DMpart4.weight(index4);

                        ComplexType MatrixElement = index2ket_iter.value()*
                                                    index2bra_iter.value()*
                                                    O3matrix.coeff(index3,index4)*
                                                    CX4matrix.coeff(index4,index1);

                        MatrixElement *= Permutation.sign;

                        if(TwoParticleGFTerm::IsRelevant(MatrixElement)){
                            term_counter++;
                            Terms.push_back(TwoParticleGFTerm(MatrixElement,beta,
                                                              E1,E2,E3,E4,
                                                              weight1,weight2,weight3,weight4,Permutation));
                        }
                    }
                    ++index2bra_iter;
                    ++index2ket_iter;
                };
            }
        };
    }
    };
    DEBUG("filling " << Terms.size() << " elements");
    //Storage->fill(Terms);
    //Terms.clear();
// DEBUG((*Storage)(0,0,0));
}

const TwoParticleGFPart::MatsubaraContainer& TwoParticleGFPart::getMatsubaraContainer(){
    return *Storage;
}

ComplexType TwoParticleGFPart::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3)
{
    // TODO: Place this variable to a wider scope?
    ComplexType MatsubaraSpacing = I*M_PI/DMpart1.getBeta();
    long MatsubaraNumberOdd1 = 2*MatsubaraNumber1 + 1;
    long MatsubaraNumberOdd2 = 2*MatsubaraNumber2 + 1;
    long MatsubaraNumberOdd3 = 2*MatsubaraNumber3 + 1;
    ComplexType Frequency1 = MatsubaraSpacing * RealType(MatsubaraNumberOdd1);
    ComplexType Frequency2 = MatsubaraSpacing * RealType(MatsubaraNumberOdd2);
    ComplexType Frequency3 = MatsubaraSpacing * RealType(MatsubaraNumberOdd3);

    ComplexType Value = 0;
    for(std::list<TwoParticleGFTerm>::const_iterator pTerm = Terms.begin(); pTerm != Terms.end(); ++pTerm)
        Value += (*pTerm)(Frequency1,Frequency2,Frequency3);

    return Value+(*Storage)(MatsubaraNumber1,MatsubaraNumber2,MatsubaraNumber3);
}
