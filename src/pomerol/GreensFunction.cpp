#include "pomerol/GreensFunction.h"

#include <cassert>
#include <stdexcept>

namespace Pomerol {

GreensFunction::GreensFunction(const StatesClassification& S, const Hamiltonian& H,
                               const AnnihilationOperator& C, const CreationOperator& CX,
                               const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), C(C), CX(CX), DM(DM)
{
}

GreensFunction::GreensFunction(const GreensFunction& GF) :
    Thermal(GF.beta), ComputableObject(GF), S(GF.S), H(GF.H), C(GF.C), CX(GF.CX), DM(GF.DM), Vanishing(GF.Vanishing)
{
    for(auto const& part: GF.parts)
        parts.emplace_back(part);
}

void GreensFunction::prepare(void)
{
    if(getStatus() >= Prepared) return;

    // Find out non-trivial blocks of C and CX.
    MonomialOperator::BlocksBimap const& CNontrivialBlocks = C.getBlockMapping();
    MonomialOperator::BlocksBimap const& CXNontrivialBlocks = CX.getBlockMapping();

    auto Citer = CNontrivialBlocks.left.begin();
    auto CXiter = CXNontrivialBlocks.right.begin();

    while(Citer != CNontrivialBlocks.left.end() && CXiter != CXNontrivialBlocks.right.end()){
        // <Cleft|C|Cright><CXleft|CX|CXright>
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->second;
        BlockNumber CXright = CXiter->first;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(Cleft == CXright && Cright == CXleft){
            // check if retained blocks are included. If not, do not push.
            if (DM.isRetained(Cleft) || DM.isRetained(Cright)) {
                parts.emplace_back(
                              C.getPartFromLeftIndex(Cleft),
                              CX.getPartFromRightIndex(CXright),
                              H.getPart(Cright), H.getPart(Cleft),
                              DM.getPart(Cright), DM.getPart(Cleft));
            }
        }

        unsigned long CleftInt = Cleft;
        unsigned long CXrightInt = CXright;

        if(CleftInt <= CXrightInt) Citer++;
        if(CleftInt >= CXrightInt) CXiter++;
    }
    if(!parts.empty())
        Vanishing = false;

    setStatus(Prepared);
}

void GreensFunction::compute()
{
    if(getStatus() >= Computed) return;
    if(getStatus() < Prepared)
        prepare();

    if(getStatus() < Computed){
        for(auto & p : parts)
            p.compute();
    }

    setStatus(Computed);
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

} // namespace Pomerol
