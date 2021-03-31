#include "pomerol/GreensFunction.h"

namespace Pomerol{

template<bool Complex>
GreensFunction<Complex>::GreensFunction(const StatesClassification<Complex>& S,
                                        const Hamiltonian<Complex>& H,
                                        const AnnihilationOperator<Complex>& C,
                                        const CreationOperator<Complex>& CX,
                                        const DensityMatrix<Complex>& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), C(C), CX(CX), DM(DM), Vanishing(true)
{
}

template<bool Complex>
GreensFunction<Complex>::GreensFunction(const GreensFunction& GF) :
    Thermal(GF.beta), ComputableObject(GF), S(GF.S), H(GF.H), C(GF.C), CX(GF.CX), DM(GF.DM), Vanishing(GF.Vanishing)
{
    for(auto iter = GF.parts.begin(); iter != GF.parts.end(); iter++)
        parts.push_back(new PartT(**iter));
}

template<bool Complex>
GreensFunction<Complex>::~GreensFunction()
{
    for(auto iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

template<bool Complex>
void GreensFunction<Complex>::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of C and CX.
    typename FieldOperator<Complex>::BlocksBimap const& CNontrivialBlocks = C.getBlockMapping();
    typename FieldOperator<Complex>::BlocksBimap const& CXNontrivialBlocks = CX.getBlockMapping();

    typename FieldOperator<Complex>::BlocksBimap::left_const_iterator Citer = CNontrivialBlocks.left.begin();
    typename FieldOperator<Complex>::BlocksBimap::right_const_iterator CXiter = CXNontrivialBlocks.right.begin();

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
                parts.push_back(new PartT(
                              (AnnihilationOperatorPart<Complex>&)C.getPartFromLeftIndex(Cleft),
                              (CreationOperatorPart<Complex>&)CX.getPartFromRightIndex(CXright),
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

template<bool Complex>
void GreensFunction<Complex>::compute()
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
unsigned short GreensFunction<Complex>::getIndex(size_t Position) const
{
    switch(Position){
        case 0: return C.getIndex();
        case 1: return CX.getIndex();
        default: assert(0);
    }
    throw std::logic_error("GreensFunction :: wrong operator");
    return C.getIndex();
}

template<bool Complex>
bool GreensFunction<Complex>::isVanishing(void) const
{
    return Vanishing;
}

} // end of namespace Pomerol
