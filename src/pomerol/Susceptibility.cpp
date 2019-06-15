#include "pomerol/Susceptibility.h"

namespace Pomerol{

Susceptibility::Susceptibility(const StatesClassification& S, const Hamiltonian& H,
                               const AnnihilationOperator& A, const CreationOperator& B,
                               const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), B(B), DM(DM), Vanishing(true)
{
}

Susceptibility::Susceptibility(const Susceptibility& Chi) :
    Thermal(Chi.beta), ComputableObject(Chi), S(Chi.S), H(Chi.H), A(Chi.A), B(Chi.B), DM(Chi.DM), Vanishing(Chi.Vanishing)
{
    for(std::list<SusceptibilityPart*>::const_iterator iter = Chi.parts.begin(); iter != Chi.parts.end(); iter++)
        parts.push_back(new SusceptibilityPart(**iter));
}

Susceptibility::~Susceptibility()
{
    for(std::list<SusceptibilityPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

void Susceptibility::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A and B.
    FieldOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();
    FieldOperator::BlocksBimap const& BNontrivialBlocks = B.getBlockMapping();

    FieldOperator::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin();
    FieldOperator::BlocksBimap::right_const_iterator Biter = BNontrivialBlocks.right.begin();

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
                parts.push_back(new SusceptibilityPart(
                              (AnnihilationOperatorPart&)A.getPartFromLeftIndex(Aleft),
                              (CreationOperatorPart&)B.getPartFromRightIndex(Bright),
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

void Susceptibility::compute()
{
    if(Status>=Computed) return;
    if(Status<Prepared) prepare();

    if(Status<Computed){
        for(std::list<SusceptibilityPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
            (*iter)->compute();
    }
    Status = Computed;
}

unsigned short Susceptibility::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return A.getIndex();
        case 1: return B.getIndex();
        default: assert(0);
    }
    throw std::logic_error("Susceptibility :: wrong operator");
    return A.getIndex();
}

bool Susceptibility::isVanishing(void) const
{
    return Vanishing;
}

} // end of namespace Pomerol
