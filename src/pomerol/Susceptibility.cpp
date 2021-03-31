#include "pomerol/Susceptibility.h"

namespace Pomerol{

template<bool Complex>
Susceptibility<Complex>::Susceptibility(const StatesClassification<Complex>& S,
                                        const Hamiltonian<Complex>& H,
                                        const QuadraticOperator<Complex>& A,
                                        const QuadraticOperator<Complex>& B,
                                        const DensityMatrix<Complex>& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), B(B), DM(DM), Vanishing(true),
    ave_A(0), ave_B(0), SubtractDisconnected(false)
{
}

template<bool Complex>
Susceptibility<Complex>::Susceptibility(const Susceptibility<Complex>& Chi) :
    Thermal(Chi.beta), ComputableObject(Chi), S(Chi.S), H(Chi.H), A(Chi.A), B(Chi.B), DM(Chi.DM),
    Vanishing(Chi.Vanishing), ave_A(Chi.ave_A), ave_B(Chi.ave_B), SubtractDisconnected(Chi.SubtractDisconnected)
{
    for(auto iter = Chi.parts.begin(); iter != Chi.parts.end(); iter++)
        parts.push_back(new PartT(**iter));
}

template<bool Complex>
Susceptibility<Complex>::~Susceptibility()
{
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

template<bool Complex>
void Susceptibility<Complex>::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A and B.
    typename FieldOperator<Complex>::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();
    typename FieldOperator<Complex>::BlocksBimap const& BNontrivialBlocks = B.getBlockMapping();

    typename FieldOperator<Complex>::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin();
    typename FieldOperator<Complex>::BlocksBimap::right_const_iterator Biter = BNontrivialBlocks.right.begin();

    while(Aiter != ANontrivialBlocks.left.end() && Biter != BNontrivialBlocks.right.end()){
        // <Aleft|A|Aright><Bleft|B|Bright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;
        BlockNumber Bleft = Biter->second;
        BlockNumber Bright = Biter->first;


        // Select a relevant 'world stripe' (sequence of blocks).
        if(Aleft == Bright && Aright == Bleft){
        //DEBUG(S.getQuantumNumbers(Aleft) << "|" << S.getQuantumNumbers(Aright) << "||" << S.getQuantumNumbers(Bleft) << "|" << S.getQuantumNumbers(Bright) );
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) || DM.isRetained(Aright) )
                parts.push_back(new PartT(
                              (QuadraticOperatorPart<Complex>&)A.getPartFromLeftIndex(Aleft),
                              (QuadraticOperatorPart<Complex>&)B.getPartFromRightIndex(Bright),
                              H.getPart(Aright), H.getPart(Aleft),
                              DM.getPart(Aright), DM.getPart(Aleft)));
        }

        unsigned long AleftInt = Aleft;
        unsigned long BrightInt = Bright;

        if(AleftInt <= BrightInt) Aiter++;
        if(AleftInt >= BrightInt) Biter++;
    }
    if (parts.size() > 0) Vanishing = false;

    Status = Prepared;
}

template<bool Complex>
void Susceptibility<Complex>::compute()
{
    if(Status>=Computed) return;
    if(Status<Prepared) prepare();

    if(Status<Computed){
        for(auto iter = parts.begin(); iter != parts.end(); iter++)
            (*iter)->compute();
    }
    Status = Computed;
}

template<bool Complex>
void Susceptibility<Complex>::subtractDisconnected()
{
    EnsembleAverage<Complex> EA_A(S, H, A, DM);
    EnsembleAverage<Complex> EA_B(S, H, B, DM);
    subtractDisconnected(EA_A, EA_B);
}

template<bool Complex>
void Susceptibility<Complex>::subtractDisconnected(ComplexType ave_A, ComplexType ave_B)
{
    SubtractDisconnected = true;
    this->ave_A = ave_A;
    this->ave_B = ave_B;
}

template<bool Complex>
void Susceptibility<Complex>::subtractDisconnected(EnsembleAverage<Complex> &EA_A,
                                                   EnsembleAverage<Complex> &EA_B)
{
    EA_A.prepare();
    EA_B.prepare();
    subtractDisconnected(EA_A.getResult(), EA_B.getResult());
}

template<bool Complex>
bool Susceptibility<Complex>::isVanishing(void) const
{
    return Vanishing;
}

template class Susceptibility<false>;
template class Susceptibility<true>;

} // end of namespace Pomerol
