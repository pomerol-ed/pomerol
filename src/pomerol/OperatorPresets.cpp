#include "pomerol/OperatorPresets.h"

namespace Pomerol {
namespace OperatorPresets {

//
// Operator N
//

N::N(ParticleIndex Nmodes):Operator(),Nmodes(Nmodes)
{
    for (ParticleIndex index=0; index<Nmodes; ++index) {
        (*this)+=n(index);
    };
};
    
std::map<FockState,ComplexType> N::actRight(const FockState &ket) const
{
    std::map<FockState,ComplexType> output;
    output[ket]=this->getMatrixElement(ket);
    return output;
}

ComplexType N::getMatrixElement(const FockState &bra, const FockState &ket) const
{
    return (bra!=ket)?0:getMatrixElement(ket);
}

ComplexType N::getMatrixElement(const FockState &ket) const
{
    return ket.count();
}

//
// Operator Sz
//

Sz::Sz(ParticleIndex Nmodes, const std::vector<ParticleIndex> & SpinUpIndices) 
    : Operator(),Nmodes(Nmodes),SpinUpIndices(SpinUpIndices)
{
    SpinDownIndices.reserve(Nmodes/2);
    for (ParticleIndex i=0; i<Nmodes; i++) { 
        if (std::find(SpinUpIndices.begin(), SpinUpIndices.end(), i)==SpinUpIndices.end()) SpinDownIndices.push_back(i);
        };
    if (SpinUpIndices.size() != SpinDownIndices.size() ) { throw ( Pomerol::Operator::exWrongLabel() ); ERROR("Sz operator requires even number of indices"); }; 
    generateTerms();
}


Sz::Sz(const std::vector<ParticleIndex> & SpinUpIndices, const std::vector<ParticleIndex> & SpinDownIndices) 
    : Operator(),Nmodes(SpinUpIndices.size() + SpinDownIndices.size()),SpinUpIndices(SpinUpIndices), SpinDownIndices(SpinDownIndices)
{
    if (SpinUpIndices.size() != SpinDownIndices.size() ) { throw ( Pomerol::Operator::exWrongLabel() ); ERROR("Sz operator requires even number of indices"); }; 
    generateTerms();
}

void Sz::generateTerms()
{
    for (ParticleIndex i=0; i<SpinUpIndices.size(); ++i) {
        (*this)+=n(SpinUpIndices[i])*0.5;
        (*this)-=n(SpinDownIndices[i])*0.5;
    }
}

ComplexType Sz::getMatrixElement(const FockState &ket) const
{
    int up_value=0, down_value=0;
    std::vector<ParticleIndex>::const_iterator it_up=SpinUpIndices.begin(), it_down=SpinDownIndices.begin();
    for (;it_up!=SpinUpIndices.end();it_up++) up_value+=ket.test(*it_up);
    for (;it_down!=SpinDownIndices.end();it_down++) down_value+=ket.test(*it_down);
    return 0.5*(up_value-down_value);
}

ComplexType Sz::getMatrixElement(const FockState &bra, const FockState &ket) const
{
    return (bra!=ket)?0:getMatrixElement(ket);
}

std::map<FockState,ComplexType> Sz::actRight(const FockState &ket) const
{
    std::map<FockState,ComplexType> output;
    output[ket]=this->getMatrixElement(ket);
    return output;
}

} // end of namespace OperatorPresets
} // end of namespace Pomerol
