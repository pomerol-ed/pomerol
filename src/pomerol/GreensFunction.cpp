#include "pomerol/GreensFunction.h"

namespace Pomerol{

GreensFunction::GreensFunction(const StatesClassification& S, const Hamiltonian& H, 
                               const AnnihilationOperator& C, const CreationOperator& CX,
                               const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), C(C), CX(CX), DM(DM), Vanishing(true)
{
}

GreensFunction::GreensFunction(const GreensFunction& GF) :
    Thermal(GF.beta), ComputableObject(GF), S(GF.S), H(GF.H), C(GF.C), CX(GF.CX), DM(GF.DM), Vanishing(GF.Vanishing)
{
    for(std::list<GreensFunctionPart*>::const_iterator iter = GF.parts.begin(); iter != GF.parts.end(); iter++)
        parts.push_back(new GreensFunctionPart(**iter));
}

GreensFunction::~GreensFunction()
{
    for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

void GreensFunction::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of C and CX.
    FieldOperator::BlocksBimap const& CNontrivialBlocks = C.getBlockMapping();
    FieldOperator::BlocksBimap const& CXNontrivialBlocks = CX.getBlockMapping();

    FieldOperator::BlocksBimap::left_const_iterator Citer = CNontrivialBlocks.left.begin();
    FieldOperator::BlocksBimap::right_const_iterator CXiter = CXNontrivialBlocks.right.begin();

    while(Citer != CNontrivialBlocks.left.end() && CXiter != CXNontrivialBlocks.right.end()){
        // <Cleft|C|Cright><CXleft|CX|CXright>
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->second;
        BlockNumber CXright = CXiter->first;


        // Select a relevant 'world stripe' (sequence of blocks).
        if(Cleft == CXright && Cright == CXleft){
        //DEBUG(S.getQuantumNumbers(Cleft) << "|" << S.getQuantumNumbers(Cright) << "||" << S.getQuantumNumbers(CXleft) << "|" << S.getQuantumNumbers(CXright) );
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Cleft) || DM.isRetained(Cright) )
                parts.push_back(new GreensFunctionPart(
                              (AnnihilationOperatorPart&)C.getPartFromLeftIndex(Cleft),
                              (CreationOperatorPart&)CX.getPartFromRightIndex(CXright),
                              H.getPart(Cright), H.getPart(Cleft),
                              DM.getPart(Cright), DM.getPart(Cleft)));
        }

        unsigned long CleftInt = Cleft;
        unsigned long CXrightInt = CXright;

        if(CleftInt <= CXrightInt) Citer++;
        if(CleftInt >= CXrightInt) CXiter++;
    }
    if (parts.size() > 0) Vanishing = false;

    Status = Prepared;
}

void GreensFunction::compute()
{
    if(Status>=Computed) return;
    if(Status<Prepared) prepare();

    if(Status<Computed){
        for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
            (*iter)->compute();
    }
    Status = Computed;
}

unsigned short GreensFunction::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C.getIndex();
        case 1: return CX.getIndex();
        default: assert(0);
    }
    throw std::logic_error("GreensFunction :: wrong operator");
    return C.getIndex();
}

bool GreensFunction::isVanishing(void) const
{
    return Vanishing;
}

} // end of namespace Pomerol
